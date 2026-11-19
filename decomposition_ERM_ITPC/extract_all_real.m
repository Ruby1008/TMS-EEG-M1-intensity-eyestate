P = hhsa_setup(); F = dir(fullfile(P.out, 'timf_spm_spTMS_*_Real_*.mat')); ids = erase({F.name}, {'timf_spm_', '.mat'});
fid = fopen(fullfile(P.out, 'extract_all_real.log'), 'w');
for k = 1:numel(ids)
    try, R = erm_extract(ids{k}); I = itpc_extract(ids{k}); fprintf(fid, '%s ok trials %s\n', ids{k}, mat2str(R.nTrials));
    catch ME, fprintf(fid, '%s FAILED %s\n', ids{k}, ME.message); end
end
fclose(fid);
