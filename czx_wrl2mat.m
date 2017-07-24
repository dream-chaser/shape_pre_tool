%% wrl
% fid = fopen('C:\Users\ChenZhixing\Desktop\000.wrl');
% vertex = fscanf(fid,'%f %f %f,',[3,36678])';
% face = fscanf(fid,'%d, %d, %d, -1,',[3,70589])';
% X = vertex(:,1);
% Y = vertex(:,2);
% Z = vertex(:,3);
% TRIV = face+1;
% shape = struct('TRIV', TRIV, 'X', X, 'Y', Y, 'Z', Z);

%% obj + path
nv = 9469;
fid = fopen('J:\BP4D-preprocessed\F001\T8\000.obj');
vertex = fscanf(fid,'v %f %f %f\n',[3,nv])';
face = fscanf(fid,'f %d/%d/%d %d/%d/%d %d/%d/%d\n',[9,inf])';
% face = fscanf(fid,'f %d %d %d\n',[3,inf])';
% path = fscanf(fid,'%f %f %f\n',[3,inf])';
X = vertex(:,1);
Y = vertex(:,2);
Z = vertex(:,3);
TRIV = [face(:,1),face(:,4),face(:,7)];
% TRIV = [face(:,1),face(:,2),face(:,3)];
shape = struct('TRIV', TRIV, 'X', X, 'Y', Y, 'Z', Z);
fclose(fid);
clear nv fid vertex face X Y Z TRIV ans;

%% draw disk
% fid = fopen('C:\Users\ChenZhixing\Desktop\test.txt');
% t=fscanf(fid,'%d %f %f\n',[3,inf])';
% fclose(fid);
% idx = int32(t(:,1));
% rho_color=zeros(9256,1);
% theta_color=zeros(9256,1);
% rho_color(idx)=t(:,2);
% theta_color(idx)=t(:,3);
% theta_color=theta_color/max(theta_color);
% rho_color=rho_color/max(rho_color);
% czx_mat2wrl(shape,rho_color,'C:\Users\ChenZhixing\Desktop\rho3.wrl');
% czx_mat2wrl(shape,theta_color,'C:\Users\ChenZhixing\Desktop\theta3.wrl');
% clear fid t idx ans;% rho_color theta_color;

%% draw path
% pi = length(shape.X)+1;
% shape.X = [shape.X;path(:,1)];
% shape.Y = [shape.Y;path(:,2)];
% shape.Z = [shape.Z;path(:,3)];
% shape.X = [shape.X;path(:,1)+0.1];
% shape.Y = [shape.Y;path(:,2)+0.1];
% shape.Z = [shape.Z;path(:,3)+0.1];
% color = zeros(length(shape.X),1);
% color(pi) = 1/6;
% for i=1:length(path)-1
%     shape.TRIV = [shape.TRIV;[pi,pi+1,pi+length(path)]];
%     pi = pi + 1;
%     color(pi) = 5/6;
% end
% color(pi) = 1/6;
% czx_mat2obj(shape,'C:\Users\ChenZhixing\Desktop\path.obj');
% % czx_mat2wrl(shape,color,'C:\Users\ChenZhixing\Desktop\path.wrl');
% clear ans pi color i;

%% curv
% G = curv.G;
% H = curv.H;
% k1 = curv.k1;
% k2 = curv.k2;

% v_G  = zeros(length(shape.X),1);
% v_H  = zeros(length(shape.X),1);
% % v_k1 = zeros(length(shape.X),1);
% % v_k2 = zeros(length(shape.X),1);
% 
% f = shape.TRIV;
% 
% for i=1:length(shape.X)
%     fprintf('%d/%d\n',i,length(shape.X));
%     nei_f = [find(f(:,1)==i);find(f(:,2)==i);find(f(:,3)==i)];
%     v_G(i) = mean(G(nei_f));
%     v_H(i) = mean(H(nei_f));
% %     v_k1(i) = mean(k1(nei_f));
% %     v_k2(i) = mean(k2(nei_f));
% end

% fid = fopen('C:\Users\ChenZhixing\Desktop\curv.csv');
% G = fscanf(fid,'%d,%f',[2,inf])';
