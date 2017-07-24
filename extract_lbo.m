function extract_lbo(srcpath, dstpath, nLBO)

fnames = dir(fullfile(srcpath, '*.mat'));
parfor i = 1 : length(fnames)
    if exist(fullfile(dstpath, fnames(i).name), 'file')
%         fprintf('%s already processed, skipping\n', fnames(i).name)
        continue;
    end
%     fprintf('Processing %s\n', fnames(i).name);
    fprintf('%s\n', fnames(i).name);
    continue;
    tmp = load(fullfile(srcpath, fnames(i).name));
    try
        [Phi, Lambda, A] = calc_lbo(tmp.shape, nLBO);
        if(~exist(dstpath, 'dir'))
            mkdir(dstpath);
        end
        parsave(fullfile(dstpath, fnames(i).name), Phi, Lambda, A);
    catch
        fprintf('Error in Processing %s\n', fnames(i).name)
    end
end
end

function parsave(fn, Phi, Lambda, A)
save(fn, 'Phi', 'Lambda', 'A', '-v7.3')
end
