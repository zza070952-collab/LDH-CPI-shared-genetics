#!/usr/bin/env python3
"""Execute a declared P2 task with pause/resource guards and verifiable completion."""
import argparse, csv, fcntl, gzip, hashlib, json, os, resource, shutil, signal, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ADMIN = ROOT / '00_admin'

def stamp():
    return datetime.now(timezone.utc).isoformat()

def digest(path):
    path = Path(path).resolve()
    stat = path.stat()
    key = str(path)
    cache_path = ADMIN / 'INPUT_SHA256_CACHE.json'
    with (ADMIN / 'INPUT_SHA256_CACHE.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
        old = cache.get(key, {})
        identity = [stat.st_size, stat.st_mtime_ns, stat.st_ino, stat.st_ctime_ns]
        if old.get('identity') == identity:
            return old['sha256']
        h = hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda: f.read(8*1024*1024), b''):
                h.update(chunk)
        cache[key] = {'identity':identity,'sha256':h.hexdigest(),'verified_at':stamp()}
        temp = cache_path.with_suffix('.json.tmp'); temp.write_text(json.dumps(cache,indent=2)); temp.replace(cache_path)
        return h.hexdigest()

def memory_available_mb():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            return int(line.split()[1]) // 1024
    raise RuntimeError('MemAvailable missing')

def disk_reserve_bytes(spec):
    policy_path=ADMIN/'DISK_POLICY_CURRENT.json'
    if policy_path.exists():
        policy=json.loads(policy_path.read_text())
        assert policy['authorization']=='USER_EXPLICIT_CANCEL_80_GIB_20260910'
        assert policy['fixed_reserve_bytes']==0
        return 0
    return int(spec.get('min_free_disk_bytes',80*1024**3))

def resolve_root_path(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path

def validate_output(item):
    path = resolve_root_path(item['path'])
    if ROOT not in path.resolve().parents:
        raise ValueError('Output escapes new project root: '+str(path))
    if not path.is_file() or path.stat().st_size < item.get('min_bytes',1):
        raise ValueError('Missing/empty output: '+str(path))
    if item.get('columns'):
        opener = gzip.open if path.suffix == '.gz' else open
        with opener(path,'rt') as f:
            header = f.readline().rstrip('\n\r').split(item.get('delimiter','\t'))
            missing = set(item['columns']) - set(header)
            if missing: raise ValueError('Missing columns '+str(missing)+' in '+str(path))
    return {'path':str(path),'bytes':path.stat().st_size,'sha256':digest(path)}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('task'); args=parser.parse_args()
    spec_path=ROOT/'16_workflow/tasks'/f'{args.task}.json'
    spec=json.loads(spec_path.read_text())
    run_dir=ROOT/'logs'/args.task; run_dir.mkdir(parents=True,exist_ok=True)
    task_lock=(run_dir/'task.lock').open('a+')
    try:
        fcntl.flock(task_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:
        print('ALREADY_RUNNING_LOCK '+args.task);return 75
    # Also protect jobs launched by the prior runner before flock was added.
    pidfile=run_dir/'pid.json'
    if pidfile.exists():
        process=Path('/proc')/str(json.loads(pidfile.read_text())['pid'])
        try:
            if (process/'cwd').resolve()==ROOT and (process/'cmdline').read_bytes().rstrip(b'\0').split(b'\0')==[str(x).encode() for x in spec['command']]:
                print('ALREADY_RUNNING_EXACT_PID '+args.task);return 75
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    done=run_dir/'.done'; exitfile=run_dir/'.exitcode'; logfile=run_dir/(args.task+'.log')
    if (ADMIN/'PAUSE_REQUESTED').exists():
        (run_dir/'wrapper.exitcode').write_text('75\n')
        print('PAUSED_BY_USER_QUOTA_RESET'); return 75
    input_rows=[]
    for item in spec.get('inputs',[]):
        path=resolve_root_path(item)
        if not path.is_file(): raise FileNotFoundError(path)
        input_rows.append({'path':str(path),'sha256':digest(path),'bytes':path.stat().st_size})
    software=[]
    for name,cmd in spec.get('software_versions',{}).items():
        result=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True,timeout=60)
        if result.returncode: raise RuntimeError('Software check failed: '+name+' '+result.stdout+' '+result.stderr)
        software.append({'software':name,'version':(result.stdout+result.stderr).strip()})
    signature=hashlib.sha256(json.dumps({'spec':spec,'inputs':input_rows,'software':software},sort_keys=True).encode()).hexdigest()
    if done.exists():
        prior=json.loads(done.read_text())
        if prior.get('signature')==signature:
            current=[validate_output(x) for x in spec.get('outputs',[])]
            if current==prior.get('outputs'):
                (run_dir/'wrapper.exitcode').write_text('0\n')
                print('REUSE_EXACT '+args.task); return 0
        raise RuntimeError('Existing completion differs; explicit new task id/checkpoint required')
    reserve=disk_reserve_bytes(spec)
    if shutil.disk_usage(ROOT).free < reserve + int(spec.get('estimated_output_bytes',0)):
        (run_dir/'wrapper.exitcode').write_text('74\n')
        (run_dir/'BLOCKED.json').write_text(json.dumps({'reason':'DISK_RESERVE','time':stamp()},indent=2)); return 74
    if memory_available_mb()<61440+int(spec.get('mem_mb',4096)):
        (run_dir/'wrapper.exitcode').write_text('74\n')
        (run_dir/'BLOCKED.json').write_text(json.dumps({'reason':'RAM_RESERVE','time':stamp()},indent=2)); return 74
    for filename,rows,fields in [('input_checksums.tsv',input_rows,['path','sha256','bytes']),('software_versions.tsv',software,['software','version'])]:
        with (run_dir/filename).open('w',newline='') as f:
            w=csv.DictWriter(f,fields,delimiter='\t'); w.writeheader();w.writerows(rows)
    env=os.environ.copy();threads=str(spec.get('threads',1))
    env.update({'OMP_NUM_THREADS':threads,'OPENBLAS_NUM_THREADS':threads,'MKL_NUM_THREADS':threads,'NUMEXPR_NUM_THREADS':threads})
    env.update(spec.get('env',{}))
    def limits():
        cap=int(spec.get('mem_mb',4096))*1024**2
        resource.setrlimit(resource.RLIMIT_AS,(cap,cap))
    with logfile.open('a') as log:
        log.write(stamp()+' START '+json.dumps(spec['command'])+'\n');log.flush()
        proc=subprocess.Popen(spec['command'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True,preexec_fn=limits)
        (run_dir/'pid.json').write_text(json.dumps({'pid':proc.pid,'started':stamp(),'task':args.task}))
        disk_paused=False
        while proc.poll() is None:
            reserve=disk_reserve_bytes(spec)
            free=shutil.disk_usage(ROOT).free
            if free<reserve and not disk_paused:
                os.killpg(proc.pid,signal.SIGSTOP);disk_paused=True
                (run_dir/'RESOURCE_PAUSED.json').write_text(json.dumps({'reason':'DISK_RESERVE','free_bytes':free,'pid':proc.pid,'time':stamp()}))
            time.sleep(5)
        code=proc.returncode;exitfile.write_text(str(code)+'\n')
        log.write(stamp()+' EXIT '+str(code)+'\n')
    if code:
        (run_dir/'wrapper.exitcode').write_text(str(code)+'\n')
        return code
    try:
        outputs=[validate_output(x) for x in spec.get('outputs',[])]
    except Exception as e:
        (run_dir/'wrapper.exitcode').write_text('65\n')
        (run_dir/'OUTPUT_QC_FAILED.json').write_text(json.dumps({'reason':str(e),'time':stamp()}));return 65
    done.write_text(json.dumps({'task':args.task,'status':'SUCCESS','signature':signature,'actual_exit_code':code,'outputs':outputs,'completed':stamp()},indent=2))
    (run_dir/'wrapper.exitcode').write_text('0\n')
    return 0

if __name__=='__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        if len(sys.argv)>1 and sys.argv[1].replace('_','').replace('-','').isalnum():
            failure_dir=ROOT/'logs'/sys.argv[1]
            failure_dir.mkdir(parents=True,exist_ok=True)
            (failure_dir/'wrapper.exitcode').write_text('70\n')
            (failure_dir/('preflight_failure_'+str(time.time_ns())+'.json')).write_text(json.dumps({'time':stamp(),'exception':type(exc).__name__,'reason':str(exc),'algorithm_exit':'NOT_RUN_OR_SEPARATELY_RECORDED'},indent=2))
        raise
