a1 = 0/180*pi;  %À≥ ±’Î
a2 = 0/180*pi;  %”“
a3 = 0/180*pi; %…œ

R1=[cos(a1),-sin(a1),0;sin(a1),cos(a1),0;0,0,1];
R2=[cos(a2),0,-sin(a2);0,1,0;sin(a2),0,cos(a2)];
R3=[1,0,0;0,cos(a3),-sin(a3);0,sin(a3),cos(a3)];

vertex = [shape.X,shape.Y,shape.Z];
vertex = vertex*R1*R2*R3;
TRIV = shape.TRIV;
X = vertex(:,1);
Y = vertex(:,2);
Z = vertex(:,3);
newshape = struct('TRIV',TRIV,'X',X,'Y',Y,'Z',Z);
czx_mat2obj(newshape, 'C:\Users\ChenZhixing\Desktop\test.obj');

clear a1 a2 a3 R1 R2 R3 T vertex TRIV X Y Z newshape;
clc