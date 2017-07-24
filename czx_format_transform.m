function czx_format_transform()

addpath(genpath('isc'));
rad          = 0.01;    % disk radius

srcpath = 'E:\Research\Code\ShapeNet\ShapeNet_modified_by_czx\data\BU3D\shapes\test';
%datapath = 'E:\Research\datasets\BU3D\BU3DFETri_from_zqk';
dstpath = 'E:\Research\Code\ShapeNet\ShapeNet_modified_by_czx\data\BU3D\dists-error';

fnames = dir(fullfile(srcpath, '*.mat'));
parfor f_i = 1 : length(fnames)
    try
        fprintf('%s\n', fnames(f_i).name);
        tmp = load(fullfile(srcpath, fnames(f_i).name));
        shape = tmp.shape;

        idxs_  = compute_indices_FPS(50, shape, 'geod');
        dists_ = compute_geodesic_dist_1vsAll(shape,idxs_,1e+07);
        diam = max(dists_(:));
        if diam > 1e+06
            throw(MException('Id:id', 'Too large diam!'));
        end
        shape = scale_shape(shape,1/diam);

        dists = zeros(size(shape.X,1),size(shape.X,1));
        sz = size(shape.X,1);
        for i = 1:sz
%             if mod(i,100) == 0
%                 fprintf('\n    %d/%d', i, sz);
%             end
            dists(i,:) = compute_geodesic_dist_1vsAll(shape,i,rad*1.5);
        end
        shape.dists = dists;

        parsave_shape(fullfile(srcpath, fnames(f_i).name), shape);
    catch ME
        fprintf('%s\n', ME.message);
        fprintf('error: %s\n', fnames(f_i).name);
        parsave_error(fullfile(dstpath, fnames(f_i).name));
    end
end
end

function parsave_shape(fn, shape)
save(fn, 'shape', '-append');
end

function parsave_error(fn)
err = true;
save(fn, 'err', '-v7.3');
end

% save_path = '../3DLandmark_by_sj/';

% if(~exist(save_path, 'dir'))
%     mkdir(save_path);
% end
% save(fullfile(save_path, '028_shape.mat'), 'shape');
