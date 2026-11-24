% Inserted after each original prune iteration. Only persistence; no RNG draws.
if mod(k,50)==0
    checkpoint_dir=getenv('P2_PRUNE_CHECKPOINT_DIR');
    if ~isempty(checkpoint_dir)
        checkpoint_rng=rng();
        prune_chunk=prunemask(:,k-49:k);
        checkpoint_file=fullfile(checkpoint_dir,sprintf('prune_%03d.mat',k));
        save('-v7',[checkpoint_file '.tmp'],'prune_chunk','checkpoint_rng','k');
        movefile([checkpoint_file '.tmp'],checkpoint_file);
        clear prune_chunk;
    end
end
