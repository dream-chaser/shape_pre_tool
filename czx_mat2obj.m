function czx_mat2obj(shape, outfilename)
    f = fopen(outfilename, 'w');
    for i=1:length(shape.X)
        fprintf(f, 'v %f %f %f\n', shape.X(i), shape.Y(i), shape.Z(i));
    end
    for i=1:length(shape.TRIV)
        fprintf(f, 'f %d/%d/%d %d/%d/%d %d/%d/%d\n', shape.TRIV(i,1), shape.TRIV(i,1), shape.TRIV(i,1), ...
                                                     shape.TRIV(i,2), shape.TRIV(i,2), shape.TRIV(i,2), ...
                                                     shape.TRIV(i,3), shape.TRIV(i,3), shape.TRIV(i,3));
    end
    fclose(f);
end