function [M] = compute_extraction_ACNN(shape, params)

nbinst    = params.nbinst; 
nbinsth   = params.nbinsth;
nLBO      = params.nLBO;

n_t_th    = nbinst * nbinsth;
tt        = 2.^[1:nbinst];
th        = [0:nbinsth-1]/nbinsth*2*pi;

I = [];
J = [];
V = [];
start_time = tic;
A = calcVoronoiRegsCircCent(shape.TRIV, [shape.X, shape.Y, shape.Z]);
A = abs(A);
% areascaling = full(diag(A)); 
% save('areascaling.mat', 'areascaling', '-v7.3');
for thi = 1:nbinsth
    W = calcLB_ACNN(shape,th(thi));
    try
        [evecs,evals] = eigs(W, A, nLBO, 'SM');
    catch
        [evecs,evals] = eigs(W, A, nLBO, -1e-5);
    end
    evals = abs(diag(evals));
    for ti = 1:nbinst
        fprintf('    (%d,%d)/(%d,%d) %2.0fs\n',thi,ti,nbinsth,nbinst,toc(start_time));
        hks_m = abs(evecs * full(diag(exp(-tt(ti)*evals))) * evecs');
%         G = hks_m(1,:)';
%         save(sprintf('%s_%d_%d%s', 'hks_m',thi,ti,'.mat'), 'G', '-v7.3');
%         hks_m = hks_m .* repmat(areascaling',[length(areascaling),1]);
        for i = 1:size(shape.X,1)
            keep = find(hks_m(i,:)>0.99*max(hks_m(i,:)));
            if isempty(keep),[~,keep] = max(hks_m(i,:));end
            sum_hks_mi  = max(sum(hks_m(i,keep)),eps);
            row_index = n_t_th*(i-1) + nbinsth*(ti-1) + thi;
            I = [I;row_index*ones(length(keep),1)];
            J = [J;keep'];
            V = [V;hks_m(i,keep)'./sum_hks_mi];
        end
    end
end
M = sparse(I,J,V,size(shape.X,1)*n_t_th,size(shape.X,1));