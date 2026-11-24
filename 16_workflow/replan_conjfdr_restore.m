% Restore contiguous chunks plus exact RNG state. A gap is a hard error.
checkpoint_dir=getenv('P2_PRUNE_CHECKPOINT_DIR');
first_iteration=1;
if ~isempty(checkpoint_dir)
    for checkpoint_k=50:50:niter
        checkpoint_file=fullfile(checkpoint_dir,sprintf('prune_%03d.mat',checkpoint_k));
        if exist(checkpoint_file,'file')
            assert(first_iteration==checkpoint_k-49,'Checkpoint gap');
            checkpoint=load(checkpoint_file);
            assert(checkpoint.k==checkpoint_k);
            assert(isequal(size(checkpoint.prune_chunk),[nsnp 50]));
            prunemask(:,checkpoint_k-49:checkpoint_k)=checkpoint.prune_chunk;
            rng(checkpoint.checkpoint_rng);
            first_iteration=checkpoint_k+1;
            clear checkpoint;
        end
    end
end
