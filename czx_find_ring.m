function [ ring ] = czx_find_ring( shape, rn )

    face = shape.TRIV;
    nf = size(face,1);
    nv = size(shape.X,1);
    ring = cell(nv, rn);
    for i=1:nf
        if ~any(ring{face(i,1),1}==face(i,2))
            ring{face(i,1),1} = [ring{face(i,1),1},face(i,2)];
        end
        if ~any(ring{face(i,1),1}==face(i,3))
            ring{face(i,1),1} = [ring{face(i,1),1},face(i,3)];
        end
        if ~any(ring{face(i,2),1}==face(i,1))
            ring{face(i,2),1} = [ring{face(i,2),1},face(i,1)];
        end
        if ~any(ring{face(i,2),1}==face(i,3))
            ring{face(i,2),1} = [ring{face(i,2),1},face(i,3)];
        end
        if ~any(ring{face(i,3),1}==face(i,1))
            ring{face(i,3),1} = [ring{face(i,3),1},face(i,1)];
        end
        if ~any(ring{face(i,3),1}==face(i,2))
            ring{face(i,3),1} = [ring{face(i,3),1},face(i,2)];
        end
    end
    
    for i=1:nv
        flag = zeros(nv,1);
        flag(i) = 1;
        for j=1:length(ring{i,1})
            flag(ring{i,1}(j)) = 1;
        end
        for j=2:rn
            ring_pre = ring{i,j-1};
            for k=1:length(ring_pre)
                k_ring1 = ring{ring_pre(k),1};
                for q=1:length(k_ring1)
                    if flag(k_ring1(q)) == 0
                        ring{i,j} = [ring{i,j}, k_ring1(q)];
                        flag(k_ring1(q)) = 1;
                    end
                end
            end
        end
        clear flag;
    end

end

