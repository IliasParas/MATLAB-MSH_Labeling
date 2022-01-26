clear
clc
tic
% A radius around triangles
r   = 6;
r_c = 28;
offset = [0, -8, -8];
fprintf('Radius of %.1f mm minus a radius %.1f with offset (%i, %i, %i).\n', r, r_c, offset(1), offset(2), offset(3))
th1 = pi;
transform = [cos(th1), -sin(th1), 0;
             sin(th1),  cos(th1), 0;
                0,         0,     1];

gmesh_file_name = 'brain_model_3d-Ilias.msh';
gmesh = fopen(gmesh_file_name);   % The format explaind here: http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php
DATA = dlmread(gmesh_file_name,'',5,1);

%% Writes the new file
fid2 = fopen('Brain_Model.msh','w');
for i=1:5
    old = num2str(fgetl(gmesh));
    fprintf(fid2,"%s\n",old);
end
dim_POS = str2double(old);
fprintf('The number of Nodes are: %i\n', dim_POS)

Nodes = DATA(1:dim_POS,1:3)*transform;
for i=1:dim_POS
    fgetl(gmesh);
    fprintf(fid2,"%i %f %f %f\n", i, Nodes(i,1), Nodes(i,2), Nodes(i,3));
end

for i=1:3
    old = num2str(fgetl(gmesh));   %%%%%%%%%%
    fprintf(fid2,"%s\n",old);
end
dim_TETS = str2double(old);

triangles = 0; Triangles(dim_TETS,3) = 0;
p = 1 + (dim_POS+3); i = 1; j = 0;
while DATA(p,1)==2
    triangles = triangles+1;	% counts the triangles
    fprintf(fid2,"%i 2 2 2000 2 %i %i %i\n", i, DATA(p,5), DATA(p,6), DATA(p,7));
    i = i + 1;
    p = p + 1;
    if DATA(p,3)==2001
        j = j + 1;
        Triangles(j,1:3) = DATA(p,5:7);
    end
end
Triangles = Triangles(1:j,:);
Tets      = DATA(1+dim_POS+3+triangles:end-1,5:8);

fprintf('The number of Triangles are: %i\n', triangles)
fprintf('The number of Elements are: %i\n', dim_TETS)
fprintf('Takes %0.1f sec to read the files.\n', toc) % Wait about 4 sec.

%% White matter processing
tic
x_tri = 0.33*(Nodes(Triangles(:,1),1)+Nodes(Triangles(:,2),1)+Nodes(Triangles(:,3),1));
y_tri = 0.33*(Nodes(Triangles(:,1),2)+Nodes(Triangles(:,2),2)+Nodes(Triangles(:,3),2));
z_tri = 0.33*(Nodes(Triangles(:,1),3)+Nodes(Triangles(:,2),3)+Nodes(Triangles(:,3),3));
for i=1:length(Tets)
    x_tet = 0.25*(Nodes(Tets(i,1),1)+Nodes(Tets(i,2),1)+Nodes(Tets(i,3),1)+Nodes(Tets(i,4),1));
    y_tet = 0.25*(Nodes(Tets(i,1),2)+Nodes(Tets(i,2),2)+Nodes(Tets(i,3),2)+Nodes(Tets(i,4),2));
    z_tet = 0.25*(Nodes(Tets(i,1),3)+Nodes(Tets(i,2),3)+Nodes(Tets(i,3),3)+Nodes(Tets(i,4),3));

%% Checking
    if (sum((x_tri-x_tet).^2+(y_tri-y_tet).^2+(z_tri-z_tet).^2<r^2) ~= 0) && (((x_tet-offset(1))^2+(y_tet-offset(2))^2+(z_tet-offset(3))^2 > r_c^2) ~= 0) % Only the outside triangles.
        fprintf(fid2,"%i 4 2 3002 3 %i %i %i %i\n", i+triangles, Tets(i,1), Tets(i,2), Tets(i,3), Tets(i,4));
    else
        fprintf(fid2,"%i 4 2 3001 4 %i %i %i %i\n", i+triangles, Tets(i,1), Tets(i,2), Tets(i,3), Tets(i,4));
    end
end
fprintf('Takes %0.1f sec to write the files.\n', toc) % Wait about 23 sec.

%% Finish
fprintf(fid2,"$EndElements"); % The last line.
fclose(gmesh);
fclose(fid2);
clear p x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4

