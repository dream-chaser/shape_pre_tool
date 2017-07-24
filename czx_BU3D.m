rng(42,'twister')
addpath(genpath('isc'))

windows_prefix = 'E:\Research\datasets\BU3D\BU3DFETri_from_zqk\';
linux_prefix = '~/research/';
prefix = windows_prefix;
prefix_BU3D = 'F:\BU3D\BU3D_XCNN\';
remeshed_BU3D = 'F:\BU3D\BU3D_remeshed\';
submesh_BU3D = 'D:\BU3D-submesh\';

data_path = sprintf('%s', prefix);
% ldm_path = sprintf('%s%s', data_path, 'BU3DLandmarks/');
shapes_path = sprintf('%s%s', prefix_BU3D, 'shapes');
descs_path = sprintf('%s%s', prefix_BU3D, 'descs');
objs_path = sprintf('%s%s', prefix_BU3D, 'objs');
curvs_path = sprintf('%s%s', prefix_BU3D, 'curvs');

new_shapes_path = sprintf('%s%s', remeshed_BU3D, 'shapes');

submesh_objs_path = sprintf('%s%s', submesh_BU3D, 'objs');
submesh_shapes_path = sprintf('%s%s', submesh_BU3D, 'shapes');
submesh_LBOs_path = sprintf('%s%s', submesh_BU3D, 'LBOs');
submesh_descs_path = sprintf('%s%s', submesh_BU3D, 'descs');

%    fnames = dir(fullfile(new_shapes_path, '*.mat'));
%    for i = 1 : length(fnames)
%        tmp = load(fullfile(new_shapes_path,fnames(i).name));
%        ldm_index = tmp.ldm_index;
%        f = fopen(fullfile(new_shapes_path,sprintf('%s%s', fnames(i).name(1:34), 'txt')),'w');
%        for i =1:length(ldm_index)
%             fprintf(f,'%d ',ldm_index(i));
%        end
%        fclose(f);
%    end

%% cal new landmarks
%    fnames = dir(fullfile(new_shapes_path, '*.mat'));
%    for i = 1 : length(fnames)
%        if ~exist(fullfile(shapes_path,fnames(i).name), 'file')
%            fprintf('%s\n',fnames(i).name);
%            continue;
%        end
%        if mod(i,100)==0
%            fprintf('%d/%d\n',i,length(fnames));
%        end
%        tmp = load(fullfile(shapes_path,fnames(i).name));
%        ldm = tmp.ldm_index;
%        ldm_v = [tmp.shape.X(ldm+1),tmp.shape.Y(ldm+1),tmp.shape.Z(ldm+1)];
%        tmp = load(fullfile(new_shapes_path,fnames(i).name));
%        v = [tmp.shape.X,tmp.shape.Y,tmp.shape.Z];
%        min_dist = zeros(1,length(ldm_v));
%        min_dist(:,:) = 1000;
%        ldm_index = zeros(1,length(ldm_v));
%        for j = 1 : length(v)
%            for k = 1:length(ldm_v)
%                dist = (v(j,1)-ldm_v(k,1))*(v(j,1)-ldm_v(k,1)) + ...
%                       (v(j,2)-ldm_v(k,2))*(v(j,2)-ldm_v(k,2)) + ...
%                       (v(j,3)-ldm_v(k,3))*(v(j,3)-ldm_v(k,3));
%                if dist < min_dist(1,k)
%                    min_dist(1,k)=dist;
%                    ldm_index(1,k)=j-1;
%                end
%            end
%        end
%        save(fullfile(new_shapes_path,fnames(i).name),'ldm_index','-append');
%    end


%    fnames = dir(fullfile(curvs_path, '*.mat'));
%    f_num = length(fnames);
%    for i = 1 : f_num
%        tmp = load(fullfile(curvs_path, fnames(i).name));
%        if length(tmp.shape.TRIV) ~= length(tmp.shape.curv.k1)
%            fprintf('%s\n',fnames(i).name);
%        end
%    end


%    fnames = dir(fullfile(data_path, '*.mat'));
%    f_num = length(fnames);
%    for i = 1 : f_num
%        try
%            if exist(fullfile(curvs_path, fnames(i).name), 'file')
%                 fprintf('%s already processed, skipping\n', fnames(i).name)
%                 continue
%            end
%            fprintf('%d/%d\n',i,f_num);
%            tmp = load(fullfile(data_path, fnames(i).name));
%            TRIV = tmp.face;
%            X = tmp.vertex(:,1);
%            Y = tmp.vertex(:,2);
%            Z = tmp.vertex(:,3);
%            fid = fopen(fullfile(objs_path, sprintf('%s%s', fnames(i).name(1:34), 'txt')),'r');
%            curv_data = fscanf(fid,'%d %d %d %d %f %f %f %f %f %f %f %f',[12,inf])';
%            fclose(fid);
%            curv_data = flipud(curv_data);
%            k1 = curv_data(:,5);
%            k2 = curv_data(:,6);
%            k1_dir = curv_data(:,7:9);
%            k2_dir = curv_data(:,10:12);
%            curv = struct('k1', k1, 'k2', k2, 'k1_dir', k1_dir, 'k2_dir', k2_dir);
%            shape = struct('TRIV', TRIV, 'X', X, 'Y', Y, 'Z', Z, 'curv', curv);
%            if(~exist(curvs_path, 'dir'))
%                mkdir(curvs_path);
%            end
%            parsave(fullfile(curvs_path, fnames(i).name), shape);
%        catch
%            fprintf('error: %s\n', fnames(i).name);
%        end
%    end
% 
%    function parsave(fn, shape)
%        save(fn, 'shape', '-v7.3');
%    end

% error_num = 0;
% error(1) = -1;
% error_cnt(1) = -1;
% error_file(1,:)='F0001_AN01WH_F3D_manifold_tiangle.mat';
% error_file_num = 0;
% % Female 1-56
% for i=1:56
%     file_path = sprintf('%sF%04d/', ldm_path, i);
%     shape_file_path = sprintf('%s%03d/', data_path, i);
%     %fprintf('Processing %s\n', file_path);
%     fnames = dir(fullfile(file_path, '*.mat'));
%     for j = 1 : length(fnames)
%         tmp_shape = load(fullfile(shape_file_path, [fnames(j).name(1:16),'_manifold_tiangle.mat']));
%         tmp_ldm = load(fullfile(file_path, fnames(j).name));
%         len = length(tmp_shape.vertex);
%         for k = 1:length(tmp_ldm.landmark_3D)
%             ttt = tmp_ldm.landmark_3D(k,4);
%             if ttt >= len
%                 if k <= 68
%                     file_printed = false;
%                     for ti = 1:length(error_file(:,1))
%                         if strcmp(error_file(ti,:),[fnames(j).name(1:16),'_manifold_tiangle.mat'])
%                             file_printed = true;
%                         end
%                     end
%                     if ~file_printed
%                         fprintf('%s\n', [fnames(j).name(1:16),'_manifold_tiangle.mat']);
%                         error_file_num = error_file_num + 1;
%                         error_file(error_file_num,:)=[fnames(j).name(1:16),'_manifold_tiangle.mat'];
%                     end
%                 end
%                 find_idx = find(error==k);
%                 if find_idx
%                     error_cnt(find_idx) = error_cnt(find_idx) + 1;
%                 else
%                     error_num = error_num + 1;
%                     error(error_num) = k;
%                     error_cnt(error_num) = 1;
%                 end
%             end
%         end
%     end
% end
% % Male 1-44
% for i=1:44
%     file_path = sprintf('%sM%04d/', ldm_path, i);
%     shape_file_path = sprintf('%s%03d/', data_path, i+56);
%     %fprintf('Processing %s\n', file_path);
%     fnames = dir(fullfile(file_path, '*.mat'));
%     for j = 1 : length(fnames)
%         tmp_shape = load(fullfile(shape_file_path, [fnames(j).name(1:16),'_manifold_tiangle.mat']));
%         tmp_ldm = load(fullfile(file_path, fnames(j).name));
%         len = length(tmp_shape.vertex);
%         for k = 1:length(tmp_ldm.landmark_3D)
%             ttt = tmp_ldm.landmark_3D(k,4);
%             if ttt >= len
%                 if k <= 68
%                     file_printed = false;
%                     for ti = 1:length(error_file(:,1))
%                         if strcmp(error_file(ti,:),[fnames(j).name(1:16),'_manifold_tiangle.mat'])
%                             file_printed = true;
%                         end
%                     end
%                     if ~file_printed
%                         fprintf('%s\n', [fnames(j).name(1:16),'_manifold_tiangle.mat']);
%                         error_file_num = error_file_num + 1;
%                         error_file(error_file_num,:)=[fnames(j).name(1:16),'_manifold_tiangle.mat'];
%                     end
%                 end
%                 find_idx = find(error==k);
%                 if find_idx
%                     error_cnt(find_idx) = error_cnt(find_idx) + 1;
%                 else
%                     error_num = error_num + 1;
%                     error(error_num) = k;
%                     error_cnt(error_num) = 1;
%                 end
%             end
%         end
%     end
% end
                
%% copy shape from original data
%    fnames = dir(fullfile(data_path, '*.mat'));
%    if(~exist(shapes_path, 'dir'))
%        mkdir(shapes_path);
%    end
%    for i = 1 : length(fnames)
%        if ~exist(fullfile(shapes_path,fnames(i).name), 'file')
%            continue;
%        end
%        fprintf('%d/%d\n',i,length(fnames));
%        tmp = load(fullfile(data_path, fnames(i).name));
%        TRIV = tmp.face;
%        X = tmp.vertex(:,1);
%        Y = tmp.vertex(:,2);
%        Z = tmp.vertex(:,3);
%        shape = struct('TRIV', TRIV, 'X', X, 'Y', Y, 'Z', Z);
%        tmp2 = load(fullfile(shapes_path, fnames(i).name));
%        label = tmp2.label;
%        ldm_index = tmp2.ldm_index;
%        save(fullfile(shapes_path, fnames(i).name), 'shape', 'label', 'ldm_index', '-v7.3');
%    end

%% obj to mat & label
%    fnames = dir(fullfile(submesh_objs_path, '*.obj'));
%    if(~exist(submesh_shapes_path, 'dir'))
%        mkdir(submesh_shapes_path);
%    end
%    dot_pos = 37;
%    for i = 1 : length(fnames)
%        fprintf('%d/%d\n',i,length(fnames));
% %        class_name = fnames(i).name(7:8);
% %         if class_name == 'AN'
% %             label = [1,0,0,0,0,0,0];
% %         elseif class_name == 'DI'
% %             label = [0,1,0,0,0,0,0];
% %         elseif class_name == 'FE'
% %             label = [0,0,1,0,0,0,0];
% %         elseif class_name == 'HA'
% %             label = [0,0,0,1,0,0,0];
% %         elseif class_name == 'NE'
% %             label = [0,0,0,0,1,0,0];
% %         elseif class_name == 'SA'
% %             label = [0,0,0,0,0,1,0];
% %         elseif class_name == 'SU'
% %             label = [0,0,0,0,0,0,1];
% %         end
%         shape = czx_obj__read(fullfile(submesh_objs_path,fnames(i).name));
%         
%         save(fullfile(submesh_shapes_path, sprintf('%s%s',fnames(i).name(1:dot_pos),'mat')), 'shape', '-v7.3');
%    end


%% mat to obj
%    fnames = dir(fullfile(shapes_path, '*.mat'));
%    if(~exist(objs_path, 'dir'))
%        mkdir(objs_path);
%    end
%    for i = 1 : length(fnames)
%        fprintf('%d/%d\n',i,length(fnames));
%        tmp = load(fullfile(shapes_path, fnames(i).name));
%        shape = tmp.shape;
%        czx_mat2obj(shape,fullfile(objs_path, sprintf('%s%s',fnames(i).name(1:34),'obj')));
%    end

% % Compute LBO
% nLBO = 30;
% if(~exist(submesh_LBOs_path, 'dir'))
%     mkdir(submesh_LBOs_path);
% end
% extract_lbo(submesh_shapes_path, submesh_LBOs_path, nLBO);

% % Compute GEOVEC
% nGEOVEC = 30;
% if(~exist(submesh_descs_path, 'dir'))
%     mkdir(submesh_descs_path);
% end
% geovec_params = estimate_geovec_params(submesh_LBOs_path, nGEOVEC);
% %save('geovec_params.mat', 'geovec_params');
% extract_geovec(submesh_LBOs_path, submesh_descs_path, geovec_params);

%% label and landmarks
%     fnames = dir(fullfile(shapes_path, '*.mat'));
%     for j = 1 : length(fnames)
%         tmp = load(fullfile('F:\\BU3D_XCNN\\shapes\\', fnames(j).name));
%         ldm_index = tmp.ldm_index;
%         label = tmp.label;
%         if(~exist(shapes_path, 'dir'))
%             mkdir(shapes_path);
%         end
%         full_name = fullfile(shapes_path, fnames(j).name);
%         if exist(full_name, 'file')
%             save(full_name, 'ldm_index', 'label', '-append');
%         end
%     end

% % Female 1-56
% for i=1:56
%     file_path = sprintf('%sF%04d/', ldm_path, i);
%     fprintf('Processing %s\n', file_path);
%     fnames = dir(fullfile(file_path, '*.mat'));
%     for j = 1 : length(fnames)
%         tmp = load(fullfile(file_path, fnames(j).name));
%         ldm_index = tmp.landmark_3D(1:68,4)'; % 只取前68个
%         
%         class_name = fnames(j).name(7:8);
%         if class_name == 'AN'
%             label = [1,0,0,0,0,0,0];
%         elseif class_name == 'DI'
%             label = [0,1,0,0,0,0,0];
%         elseif class_name == 'FE'
%             label = [0,0,1,0,0,0,0];
%         elseif class_name == 'HA'
%             label = [0,0,0,1,0,0,0];
%         elseif class_name == 'NE'
%             label = [0,0,0,0,1,0,0];
%         elseif class_name == 'SA'
%             label = [0,0,0,0,0,1,0];
%         elseif class_name == 'SU'
%             label = [0,0,0,0,0,0,1];
%         end
%         
%         if(~exist(shapes_path, 'dir'))
%             mkdir(shapes_path);
%         end
%         full_name = fullfile(shapes_path, [fnames(j).name(1:16),'_manifold_tiangle.mat']);
%         if exist(full_name, 'file')
%             tmp = load(full_name);
%             A = full(tmp.A);
%             save(full_name, 'A', 'ldm_index', 'label', '-append');
%         end
%     end
% end
% 
% % Male 1-44
% for i=1:44
%     file_path = sprintf('%sM%04d/', ldm_path, i);
%     fprintf('Processing %s\n', file_path);
%     fnames = dir(fullfile(file_path, '*.mat'));
%     for j = 1 : length(fnames)
%         tmp = load(fullfile(file_path, fnames(j).name));
%         ldm_index = tmp.landmark_3D(1:68,4)'; % 只取前68个
%         
%         class_name = fnames(j).name(7:8);
%         if class_name == 'AN'
%             label = [1,0,0,0,0,0,0];
%         elseif class_name == 'DI'
%             label = [0,1,0,0,0,0,0];
%         elseif class_name == 'FE'
%             label = [0,0,1,0,0,0,0];
%         elseif class_name == 'HA'
%             label = [0,0,0,1,0,0,0];
%         elseif class_name == 'NE'
%             label = [0,0,0,0,1,0,0];
%         elseif class_name == 'SA'
%             label = [0,0,0,0,0,1,0];
%         elseif class_name == 'SU'
%             label = [0,0,0,0,0,0,1];
%         end
%         
%         if(~exist(shapes_path, 'dir'))
%             mkdir(shapes_path);
%         end
%         full_name = fullfile(shapes_path, [fnames(j).name(1:16),'_manifold_tiangle.mat']);
%         if exist(full_name, 'file')
%             tmp = load(full_name);
%             A = full(tmp.A);
%             save(full_name, 'A', 'ldm_index', 'label', '-append');
%         end
%     end
% end