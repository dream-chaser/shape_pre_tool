function [ shape,color ] = ShowCurvDir( fname, k, scale, veps )

    tmp = load(fname);
    shape = tmp.shape;
    Vertices=[shape.X, shape.Y, shape.Z];
    face_v1 = Vertices(shape.TRIV(:, 1), :);
    face_v2 = Vertices(shape.TRIV(:, 2), :);
    face_v3 = Vertices(shape.TRIV(:, 3), :);
    face_v  = (face_v1+face_v2+face_v3)/3;
    e1 = face_v2 - face_v1;
    e2 = face_v3 - face_v1;
    normal = cross(e1,e2,2);
%     normal = cross(shape.curv.k1_dir,shape.curv.k2_dir,2);
    if k == 0
        dir = normal;
    elseif k == 1
        dir = shape.curv.k1_dir;
    elseif k == 2
        dir = shape.curv.k2_dir;
    end
    pcnt = length(shape.X);
    fcnt = length(shape.TRIV);
    X = zeros(fcnt*3,1);
    Y = zeros(fcnt*3,1);
    Z = zeros(fcnt*3,1);
    new_TRIV = zeros(fcnt,3,'uint32');
    color = zeros(pcnt+fcnt*3,1);
    pind = 0;
    for i = 1:fcnt
        pind = pind + 1;
        X(pind) = face_v(i,1);
        Y(pind) = face_v(i,2);
        Z(pind) = face_v(i,3);
        color(pcnt+pind) = 5/6;
        p1 = face_v(i,:) + scale*dir(i,:)/norm(dir(i,:));
        pind = pind + 1;
        X(pind) = p1(1);
        Y(pind) = p1(2);
        Z(pind) = p1(3);
        color(pcnt+pind) = 5/6;
        pind = pind + 1;
        X(pind) = p1(1)+veps;
        Y(pind) = p1(2);
        Z(pind) = p1(3);
        color(pcnt+pind) = 5/6;
        new_TRIV(i,1) = pcnt+pind;
        new_TRIV(i,2) = pcnt+pind-1;
        new_TRIV(i,3) = pcnt+pind-2;
    end
    shape.X = [shape.X;X];
    shape.Y = [shape.Y;Y];
    shape.Z = [shape.Z;Z];
    shape.TRIV = [shape.TRIV;new_TRIV];

end