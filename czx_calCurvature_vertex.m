function [ pd,curv ] = czx_calCurvature_vertex(shape)

    %shape = load('E:\Research\Code\ShapeNet\ShapeNet_modified_by_czx\data\BU3D\shapes\F0001_AN01WH_F3D_manifold_tiangle.mat').shape;
    triv = shape.TRIV;
    v = [shape.X, shape.Y, shape.Z];
    n_v = length(v);
    pd = cell(n_v,1);
    f_k1 = zeros(n_v,1);
    f_k2 = zeros(n_v,1);
    f_G = zeros(n_v,1);
    f_H = zeros(n_v,1);
    f_inner = zeros(n_v,1);
    for i=1:100
        try
            fprintf('%d/%d\n',i,n_v);
%             f_nei = [find(triv(:,1)==i);find(triv(:,2)==i);find(triv(:,3)==i)];
%             normal_v = [0,0,0];
%             for f_nei_i=1:length(f_nei)
%                 e1 = v(triv(f_nei(f_nei_i),1),:);
%                 e2 = v(triv(f_nei(f_nei_i),2),:);
%                 normal = cross(e1,e2);
%                 normal = normal / norm(normal);
%                 normal_v = normal_v + normal;
%             end
%             normal_v = normal_v / norm(normal_v);
%             rot_mat = cal_rot_mat(normal_v); % rotate to (0,0,1)'s direction
            neighbor = find_v_nei_ring(triv,v,i,2);
            v_nei = v(neighbor,:);
%             v_nei = v_nei * rot_mat';
            [~,B,C,D,E,F] = fit_surface(v_nei);
            v0 = v(i,:);
%             v0 = v0 * rot_mat';
            dx = B + 2*D*v0(1) + E*v0(2);
            dy = C + E*v0(1) + 2*F*v0(2);
            dxy = E;
            dxx = 2*D;
            dyy = 2*F;
            e = 1 + dx*dx;
            f = dx*dy;
            g = 1 + dy*dy;
            eg_f2 = 1 + dx*dx + dy*dy;
            sqrt_eg_f2 = sqrt(eg_f2);
            l = dxx/sqrt_eg_f2;
            m = dxy/sqrt_eg_f2;
            n = dyy/sqrt_eg_f2;
            aa = n*f-m*g;
            bb = n*e-g*l;
            cc = m*e-f*l;
            t1 = (-bb-sqrt(bb*bb-4*aa*cc))/(2*aa);
            t2 = (-bb+sqrt(bb*bb-4*aa*cc))/(2*aa);
            prin_dir = [1/sqrt(1+t1*t1), 1/sqrt(1+t2*t2); t1/sqrt(1+t1*t1), t2/sqrt(1+t2*t2)];
%             k1 = (l+2*m*t1+n*t1*t1)/(e+2*f*t1+g*t1*t1);
%             k2 = (l+2*m*t2+n*t2*t2)/(e+2*f*t2+g*t2*t2);
%             G = k1*k2;
%             H = (k1+k2)/2;
            G = (l*n-m*m)/eg_f2;
            H = (l*g-2*m*f+n*e)/(2*eg_f2);
            k1 = H + sqrt(H*H-G);
            k2 = H - sqrt(H*H-G);
            p_dir_3 = [1,0;0,1;dx,dy] * prin_dir;
%             p_dir_3 = rot_mat \ p_dir_3;
            pd{i} = p_dir_3;
            f_k1(i) = k1;
            f_k2(i) = k2;
            f_G(i) = G;
            f_H(i) = H;
            f_inner(i) = dot(p_dir_3(:,1),p_dir_3(:,2));
        catch
            fprintf('error: %d\n',i);
            pd{i} = [];
            f_k1(i) = 0;
            f_k2(i) = 0;
            f_G(i) = 0;
            f_H(i) = 0;
            f_inner(i) = 0;
        end
    end
    curv = struct('G', f_G, 'H', f_H, 'k1', f_k1, 'k2', f_k2, 'inner', f_inner);
end


function [neighbor] = find_v_nei_ring(triv,v,i,max_depth)
    depth = 0;
    next_depth_num = 0;
    flag = zeros(1,length(v)); % 1 if added to queue; 0 otherwise.
    queue = [];
    q_front = 1; q_tail = 1;
    queue(q_tail) = i; q_tail = q_tail + 1; next_depth_num = next_depth_num + 1; flag(i) = 1;
    cur_depth_num = next_depth_num;
    next_depth_num = 0; visit_num = 0;
    while q_front < q_tail
        cur_v = queue(q_front); q_front = q_front + 1; visit_num = visit_num + 1;
        if depth == max_depth
            continue
        end
        for i=1:3
            faces = find(triv(:,i)==cur_v);
            for j=1:length(faces)
                for k = 1:3
                    if k ~= i
                        tmp_v = triv(faces(j),k);
                        if flag(tmp_v) == 0
                            queue(q_tail) = tmp_v; q_tail = q_tail + 1; next_depth_num = next_depth_num + 1; flag(tmp_v) = 1;
                        end
                    end
                end
            end
        end
        if visit_num == cur_depth_num
            cur_depth_num = next_depth_num;
            next_depth_num = 0; visit_num = 0;
            depth = depth + 1;
        end
    end
    neighbor = queue;
end

function [ rot_mat ] = cal_rot_mat(normal) % rotate to (0,0,1)
    x = normal(1);
    y = normal(2);
    z = normal(3);
    rot_n = [0,0,1];
    normal_project_y = [x,0,z];
    theta_y = acos(dot(rot_n,normal_project_y) / norm(normal_project_y));
    if x < 0
        theta_y = -theta_y;
    end
    rot_mat_y = [cos(theta_y),0,-sin(theta_y);0,1,0;sin(theta_y),0,cos(theta_y)];
    normal_rot_y = [0,y,norm(normal_project_y)];
    theta_x = acos(dot(rot_n,normal_rot_y) / norm(normal_rot_y));
    if y < 0
        theta_x = -theta_x;
    end
    rot_mat_x = [1,0,0;0,cos(theta_x),-sin(theta_x);0,sin(theta_x),cos(theta_x)];
    rot_mat = rot_mat_x * rot_mat_y;
end

function [A,B,C,D,E,F] = fit_surface(v_nei)
    x  = 0;  
    y  = 0;  
    x2 = 0; x3  = 0;         x4  = 0;
    xy = 0; x2y = 0;         x3y = 0;
    y2 = 0; xy2 = 0; y3 = 0; x2y2= 0; xy3 = 0; y4  = 0;
    z  = 0; xz  = 0; yz = 0; x2z = 0; xyz = 0; y2z = 0;
    n = length(v_nei);
    for i=1:n
        vx = v_nei(i,1); vy = v_nei(i,2); vz = v_nei(i,3);
        x  = x + vx;  
        y  = y + vy;  
        x2 = x2+vx*vx; x3  = x3+vx*vx*vx;                    x4  = x4+vx*vx*vx*vx;
        xy = xy+vx*vy; x2y = x2y+vx*vx*vy;                   x3y = x3y+vx*vx*vx*vy;
        y2 = y2+vy*vy; xy2 = xy2+vx*vy*vy; y3 = y3+vy*vy*vy; x2y2= x2y2+vx*vx*vy*vy; xy3 = xy3+vx*vy*vy*vy; y4  = y4+vy*vy*vy*vy;
        z  = z+vz;     xz  = xz+vx*vz;     yz = yz+vy*vz;    x2z = x2z+vx*vx*vz;     xyz = xyz+vx*vy*vz;    y2z = y2z+vy*vy*vz;
    end
    mleft  = [n,  x,   y,   x2,   xy,   y2;
              x,  x2,  xy,  x3,   x2y,  xy2;
              y,  xy,  y2,  x2y,  xy2,  y3;
              x2, x3,  x2y, x4,   x3y,  x2y2;
              xy, x2y, xy2, x3y,  x2y2, xy3;
              y2, xy2, y3,  x2y2, xy3,  y4];
    mright = [z;  xz;  yz;  x2z;  xyz;  y2z];
    coef = mleft \ mright;
%     mleft  = zeros(length(v_nei),6);
%     mright = zeros(length(v_nei),1);
%     for i=1:length(v_nei)
%         mleft(i,1) = 1;
%         mleft(i,2) = v_nei(i,1);
%         mleft(i,3) = v_nei(i,2);
%         mleft(i,4) = v_nei(i,1) * v_nei(i,1);
%         mleft(i,5) = v_nei(i,1) * v_nei(i,2);
%         mleft(i,6) = v_nei(i,2) * v_nei(i,2);
%         mright(i,1)= v_nei(i,3);
%     end
%     coef = (mleft'*mleft) \ (mleft'*mright);
    A = coef(1);
    B = coef(2);
    C = coef(3);
    D = coef(4);
    E = coef(5);
    F = coef(6);
end