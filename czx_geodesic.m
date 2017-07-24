function [ dists, paths ] = czx_geodesic( triv, vertex, src )

    nv = length(vertex);
    dists = inf(nv,1);
    dists(src) = 0;
    state = zeros(nv,1); % 0,far;1,open;2,dead
    state(src) = 1;
    while ~isempty(find(state==0,1))
        min_d = inf;
        min_i = -1;
        for i=1:nv
            if state(i)==1 && dists(i)<min_d
                min_d = dists(i);
                min_i = i;
            end
        end
        state(min_i) = 2;
        ring1 = find_ring1(triv,min_i);
        for i=1:length(ring1)
            if state(ring1(i))==0
                state(ring1(i))=1;
            end
        end
        for i=1:length(ring1)
            current_vi = ring1(i);
            if state(current_vi)==1
                min_fd = inf;
                min_fi = -1;
                ring1_face = find_ring1_face(triv,current_vi);
                for j=1:length(ring1_face)
                    fi = ring1_face(j);
                    vi0 = current_vi;
                    vi_tmp = find(triv(fi,:)~=vi0);
                    vi1 = triv(fi,vi_tmp(1));
                    vi2 = triv(fi,vi_tmp(2));
                    v0 = vertex(vi0,:);
                    v1 = vertex(vi1,:);
                    v2 = vertex(vi2,:);
                    ux1 = dists(vi1);
                    ux2 = dists(vi2);
                    d = eikonal(ux1,ux2,v0,v1,v2);
                    if d < min_fd
                        min_fd = d;
                        min_fi = j;
                    end
                end
                dists(current_vi) = min(dists(current_vi),min_fd);
            end
        end
    end

end

function [ring1,ring1_face] = find_ring1(triv,i)
    ring1_face = find_ring1_face(triv,i);
    ring1_tmp = unique(triv(ring1_face,:));
    ring1 = ring1_tmp(ring1_tmp~=i);
end

function [ring1_face] = find_ring1_face(triv,i)
    ring1_face = [find(triv(:,1)==i);find(triv(:,2)==i);find(triv(:,3)==i)];
end

function [d] = eikonal(ux1,ux2,v0,v1,v2)
    if isinf(ux1) && isinf(ux2)
        d = inf;
    elseif isinf(ux1)
        d = ux2 + norm(v0 - v2);
    elseif isinf(ux2)
        d = ux1 + norm(v0 - v1);
    else
        A = [1;1];
        B = -[ux1;ux2];
        M = [v0-v1;v0-v2];
        Q = M*M';
        a = A'/Q*A;
        b = 2*A'/Q*B;
        c = B'/Q*B-1;
        det = b*b-4*a*c;
        if det < 0
            %d = min(ux1+norm(v0-v1),ux2+norm(v0-v2));
            d = min(ux1,ux2)+1;
        else
            d = (-b+sqrt(det))/(2*a);
        end
    end
end