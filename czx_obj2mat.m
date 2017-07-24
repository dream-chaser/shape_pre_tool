function [shape] = czx_obj2mat(fname)
    fid = fopen(fname);
    num = fscanf(fid,'# %d %d\n', [2,1])';

    vertex = fscanf(fid,'v %f %f %f\n',[3,num(1)])';
    X = vertex(:,1);
    Y = vertex(:,2);
    Z = vertex(:,3);

    % face = fscanf(fid,'f %d/%d/%d %d/%d/%d %d/%d/%d\n',[9,num(2)])';
    % TRIV = [face(:,1),face(:,4),face(:,7)];

    face = fscanf(fid,'f %d %d %d\n',[3,num(2)])';
    TRIV = [face(:,1),face(:,2),face(:,3)];

    shape = struct('TRIV', TRIV, 'X', X, 'Y', Y, 'Z', Z);
    fclose(fid);
end