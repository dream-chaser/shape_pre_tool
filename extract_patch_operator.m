function extract_patch_operator(srcpath, dstpath, patch_params)

if ~exist(dstpath, 'dir')
    mkdir(dstpath);
end

fnames = dir(fullfile(srcpath, '*.mat'));
for i = 1 : length(fnames)
%     try
        if exist(fullfile(dstpath, fnames(i).name), 'file')
            fprintf('%s already processed, skipping\n', fnames(i).name)
            continue
        end
        fprintf('Processing %s\n', fnames(i).name)
        tmp = load(fullfile(srcpath, fnames(i).name));
        shape = tmp.shape;
    
        %[M, ~] = compute_extraction(shape, patch_params);
        M = compute_extraction_ACNN(shape, patch_params);
        % make a big matrix out of all the various M_i
        % each matrix in the cell array is stacked row after row.
        % this allows a more efficient moltiplication and handling in theano
        %M = sparse(cat(1, M{:}));
    
        parsave(fullfile(dstpath, fnames(i).name), M);
%     catch
%         fprintf('error: %s\n', fnames(i).name);
%     end
end
end

function parsave(fn, M)
save(fn, 'M', '-v7.3')
end
