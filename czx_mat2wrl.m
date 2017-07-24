function czx_mat2wrl(shape, color, filename)
f = fopen(filename, 'w');
fprintf(f, 'Shape {\n');
fprintf(f, '\tappearance Appearance {\n');
fprintf(f, '\t\tmaterial Material {\n');
fprintf(f, '\t\t\tdiffuseColor 0.5882 0.5882 0.5882\n');
fprintf(f, '\t\t\tambientIntensity 1.0\n');
fprintf(f, '\t\t\tspecularColor 0 0 0\n');
fprintf(f, '\t\t\tshininess 0.145\n');
fprintf(f, '\t\t\ttransparency 0\n');
fprintf(f, '\t\t}\n');
fprintf(f, '\t}\n');
fprintf(f, '\tgeometry IndexedFaceSet {\n');
fprintf(f, '\t\tccw TRUE\n');
fprintf(f, '\t\tsolid TRUE\n');
fprintf(f, '\t\tcoord Coordinate { point [\n');
% print point
for i=1:length(shape.X)
    fprintf(f, '\t\t\t%f,%f,%f,\n', shape.X(i), shape.Y(i), shape.Z(i));
end
fprintf(f, '\t\t\t]\n');
fprintf(f, '\t\t}\n');
fprintf(f, '\t\tcolor Color { color [\n');
% print color
% G_ind = geo_25_H(:,1);
% G_curv = geo_25_H(:,2);
G_curv = color;
std_G = std(G_curv);
mean_G = mean(G_curv);
max_G = max(G_curv);
min_G = min(G_curv);
for i=1:length(shape.X)
%     [r,g,b] = deal(1.0);
%     find_num = 0;
%     if find(predict_anm_index==i)
%         [g,b] = deal(0.0);
%         find_num = find_num + 1;
%     end
%     if find(predict_lbf_index==i)
%         [r,g] = deal(0.0);
%         find_num = find_num + 1;
%     end
%     if find_num == 2
%         [r,b] = deal(0.0);
%     end
%% curv
%     if find(G_ind==i)
%         j = find(G_ind==i);
%         tG = G_curv(j);
%     else
%         tG = 0;
%     end
%     tG = G_curv(i);
%     t_color = (tG-min_G)/(max_G-min_G);
%     %t_color = (tG-mean_G+std_G)/(2*std_G);
%     if t_color < 0
%         t_color = 0;
%     elseif t_color > 1
%         t_color = 1;
%     end
    t_color = G_curv(i);
    if t_color == 0
        [r,g,b] = deal(1.0);
    else
        [r,g,b] = float2RGB(t_color);
    end
%     if t_color < 0.5
%         t_color = t_color * 2;
%         r = 1-t_color;
%         g = t_color;
%         b = 0;
%     else
%         t_color = t_color * 2 - 1;
%         r = 0;
%         g = 1-t_color;
%         b = t_color;
%     end
    fprintf(f, '\t\t\t%f,%f,%f,\n', r, g, b);
end
fprintf(f, '\t\t\t]\n');
fprintf(f, '\t\t}\n');
fprintf(f, '\t\tcolorIndex   [\n');
% print face
for i=1:length(shape.TRIV)
    fprintf(f, '\t\t\t%d,%d,%d,-1,\n', shape.TRIV(i,1)-1, shape.TRIV(i,2)-1, shape.TRIV(i,3)-1);
end
fprintf(f, '\t\t]\n');
fprintf(f, '\t\tcreaseAngle 1.5\n');
fprintf(f, '\t\tcolorPerVertex TRUE\n');
fprintf(f, '\t\tcoordIndex [\n');
% print face
for i=1:length(shape.TRIV)
    fprintf(f, '\t\t\t%d,%d,%d,-1,\n', shape.TRIV(i,1)-1, shape.TRIV(i,2)-1, shape.TRIV(i,3)-1);
end
fprintf(f, '\t\t]\n');
fprintf(f, '\t}\n');
fprintf(f, '}');

fclose(f);
end