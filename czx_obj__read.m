function [shape]= czx_obj__read( input_file_name )
    fid = fopen(input_file_name,'r');
    vn = 0;
    while 1
        line = fgetl(fid);
        if line == -1
            break;
        end
        if line(1) == 'v'
            vn = vn + 1;
        end
        if line(1) == 'f'
            break;
        end
    end
    fclose(fid);
    
    fid = fopen(input_file_name,'r');
    
    vertex = fscanf(fid,'v %f %f %f\n',[3,vn])';
    X = vertex(:,1);
    Y = vertex(:,2);
    Z = vertex(:,3);

    % face = fscanf(fid,'f %d/%d/%d %d/%d/%d %d/%d/%d\n',[9,num(2)])';
    % TRIV = [face(:,1),face(:,4),face(:,7)];

    face = fscanf(fid,'f %d %d %d\n',[3,inf])';
    TRIV = [face(:,1),face(:,2),face(:,3)];

    shape = struct('TRIV', TRIV, 'X', X, 'Y', Y, 'Z', Z);
    fclose(fid);
end