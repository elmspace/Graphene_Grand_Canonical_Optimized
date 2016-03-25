% Program IsoSurface, run by typing run("./main.m")


% Clear all previous variables
clear all;
% Getting the coordiante variables from the data, and making them
% into the correct format, mesh
A=importdata('./xyz.dat');
%getting the concentration from the data file
B=importdata('./ABCD_Junc.dat');

[X,Y,Z]=meshgrid(A(:,2),A(:,1),A(:,3));


% The size is used for the for loop
x_size=size(A(:,2),1);
y_size=size(A(:,1),1);
z_size=size(A(:,3),1);

% Dividing up the concentration values (JB1C JCB2 JB2A JAB3)

VJB1C=zeros(x_size,y_size,z_size);
VJCB2=zeros(x_size,y_size,z_size);
VJB2A=zeros(x_size,y_size,z_size);
VJAB3=zeros(x_size,y_size,z_size);


% Taking the concentration values and putting them into the correct
% format, mesh-format
ii=1;
for i=1:x_size,
    for j=1:y_size,
        for k=1:z_size,
            VJB1C(i,j,k)=B(ii,1);
            VJCB2(i,j,k)=B(ii,2);
            VJB2A(i,j,k)=B(ii,3);
            VJAB3(i,j,k)=B(ii,4);
            ii=ii+1;
        end
    end
end
% clearing useless variables
clear i;
clear j;
clear k;
clear ii;
% Thres3old for the isosurf, 
cut=0.05;

axis vis3d;
view([1,0.3,0.2]);
axis off;
daspect('mode');


% VJB1C    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pVJB1C = patch(isosurface(X,Y,Z,VJB1C,cut),'FaceColor','red','EdgeColor','none');
qVJB1C = patch(isocaps(X,Y,Z,VJB1C,cut),'FaceColor','red','EdgeColor', ...
         'none');
alpha(pVJB1C,0.5);
alpha(qVJB1C,0.5);


% VJCB2 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pVJCB2 = patch(isosurface(X,Y,Z,VJCB2,cut),'FaceColor','blue','EdgeColor','none');
qVJCB2 = patch(isocaps(X,Y,Z,VJCB2,cut),'FaceColor','blue','EdgeColor', ...
          'none');
alpha(pVJCB2,0.5);
alpha(qVJCB2,0.5);



% VJB2A    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pVJB2A = patch(isosurface(X,Y,Z,VJB2A,cut),'FaceColor','green','EdgeColor','none');
qVJB2A = patch(isocaps(X,Y,Z,VJB2A,cut),'FaceColor','green','EdgeColor', ...
           'none');
alpha(pVJB2A,0.7);
alpha(qVJB2A,0.7);


% VJAB3    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pVJAB3 = patch(isosurface(X,Y,Z,VJAB3,cut),'FaceColor','green','EdgeColor','none');
qVJAB3 = patch(isocaps(X,Y,Z,VJAB3,cut),'FaceColor','green','EdgeColor', ...
           'none');
alpha(pVJAB3,0.7);
alpha(qVJAB3,0.7);

