% rng(42,'twister')
% addpath(genpath('isc'))
% 
%% Compute LBO
nLBO = 300;
extract_lbo('F:\BU3D\BU3D_remeshed\shapes', 'F:\BU3D\BU3D_remeshed\LBOs', nLBO);

%% Compute GEOVEC
%% train and test should be divied!!!
nGEOVEC = 150;
geovec_params = estimate_geovec_params('F:\BU3D\BU3D_remeshed\LBOs', nGEOVEC);
extract_geovec('F:\BU3D\BU3D_remeshed\LBOs', 'F:\BU3D\BU3D_remeshed\descs', geovec_params);
% 
% %% Compute patch operator (disk)
% patch_params.rad          = 0.01;    % disk radius
% patch_params.flag_dist    = 'fmm';   % possible choices: 'fmm' or 'min'
% patch_params.nbinsr       = 5;       % number of rings
% patch_params.nbinsth      = 16;      % number of rays
% patch_params.fhs          = 2.0;     % factor determining hardness of scale quantization
% patch_params.fha          = 0.01;    % factor determining hardness of angle quantization
% patch_params.geod_th      = true;
% extract_patch_operator('data/train/shapes', 'data/train/disk', patch_params);
% extract_patch_operator('data/test/shapes', 'data/test/disk', patch_params);

%% ACNN
% patch_params.nbinst       = 5; 
% patch_params.nbinsth      = 16;
% patch_params.nLBO         = 100;
% extract_patch_operator('F:\BU3D_XCNN\curvs', 'C:\Users\ChenZhixing\Desktop\tmp', patch_params);


% for i=1:16
%     for j=1:5
%         tmp=load(sprintf('%s_%d_%d%s','hks_m',i,j,'.mat'));
%         G{i}(:,j)=tmp.G;
%         tG=G{i}(:,j);
%         keep{i}{j}=find(tG>0.99*max(tG));
%         %czx_mat2wrl(shape, G{i}(:,j), sprintf('%s_%d_%d%s','C:\Users\ChenZhixing\Desktop\hks_m',i,j,'.wrl'));
%     end
% end
% 

% test = M(8001:8080,:)~=0;
% color = zeros(5851,1);
% for i=1:16
%     for j=1:1
%         row_ind = (j-1)*16+i;
%         t_color = (i-1)*5+j;
%         color(test(row_ind,:)',1)=t_color/80;
%     end
% end
% czx_mat2wrl(shape, color, 'C:\Users\ChenZhixing\Desktop\test.wrl');