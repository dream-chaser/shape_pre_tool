function geovec_params = estimate_geovec_params(base_path, nGEOVEC)
% base_path : path to LBO
fnames = dir(fullfile(base_path, '*.mat'));
evals = cell(length(fnames), 1);
tic
for i = 1 : length(fnames)
    if mod(i, 10) == 0
        fprintf('%f%%,', 100*i/length(fnames))
        toc
    end
    tmp = load(fullfile(base_path, fnames(i).name));
    evals{i} = tmp.Lambda;
end
evals = cat(1, evals{:});

evalRange = [0, prctile(evals,99)];
dEval = diff(evalRange) / (nGEOVEC - 1);
evalSamples = evalRange(1) : dEval : evalRange(2);

geovec_params.evalSamples = evalSamples;
geovec_params.dEval = dEval;
geovec_params.doNormalize = true;