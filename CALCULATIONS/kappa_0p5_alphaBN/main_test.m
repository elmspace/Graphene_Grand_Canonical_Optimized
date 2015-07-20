% Program IsoSurface, run by typing run("./main.m")
% Clear all previous variables
clear all;

FileName = 0.0;
File_String = num2str(FileName,8);


xyz = fullfile(['xyz_0.617565_.dat']);
ABCD = fullfile(['ABCD_0.617565_.dat']);



% Getting the coordiante variables from the data, and making them
% into the correct format, mesh
A=importdata(xyz);
%getting the concentration from the data file
B=importdata(ABCD);

[X,Y,Z]=meshgrid(A(:,2),A(:,1),A(:,3));


% The size is used for the for loop
x_size=size(A(:,1),1);
y_size=size(A(:,2),1);
z_size=size(A(:,3),1);
% Dividing up the concentration values       format from the SCFT
% program: A, C, B1, B2, B3, B4  (B1-A-B2-C-B3  + B4)
VA=zeros(x_size,y_size,z_size);
VC=zeros(x_size,y_size,z_size);
VB1=zeros(x_size,y_size,z_size);
VB2=zeros(x_size,y_size,z_size);
VB3=zeros(x_size,y_size,z_size);
VB4=zeros(x_size,y_size,z_size);
% Taking the concentration values and putting them into the correct
% format, mesh-format
Hom_Con = B(1,2);
Cop_Con = B(1,1);
ii=2;
for i=1:x_size,
    for j=1:y_size,
        for k=1:z_size,
            VA(i,j,k)=B(ii,1);
            VC(i,j,k)=B(ii,2);
            VB1(i,j,k)=B(ii,3);
            VB2(i,j,k)=B(ii,4);
            VB3(i,j,k)=B(ii,5);
            VB4(i,j,k)=B(ii,6);
            ii=ii+1;
        end
    end
end

del_A=1.5;
del_C=1.5;
del_B2=2.7;

del_B1=10.35;
del_B3=10.35;

del_B4=1.0;

pA=0.1*Cop_Con*del_A;
pC=0.1*Cop_Con*del_C;
pB1=0.355*Cop_Con*del_B1;
pB2=0.09*Cop_Con*del_B2;
pB3=0.355*Cop_Con*del_B3;

pB4=Hom_Con*del_B4;


% clearing useless variables
clear i;
clear j;
clear k;
clear ii;




%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Thres3old for the isosurf, 
cutA=pA;
cutC=pC;
cutB2=pB2; %center block

cutB1=pB1;
cutB3=pB3;

cutB4=pB4; %Homopolymer 1


axis vis3d;
view([1,0.9,0.2]);
axis off;
daspect('mode');

% B1
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pB1 = patch(isosurface(X,Y,Z,VB1,cutB1),'FaceColor','green','EdgeColor','none');
qB1 = patch(isocaps(X,Y,Z,VB1,cutB1),'FaceColor','green','EdgeColor', ...
          'none');
alpha(pB1,0.25);
alpha(qB1,0.25);

% B2
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pB2 = patch(isosurface(X,Y,Z,VB2,cutB2),'FaceColor','green','EdgeColor','none');
qB2 = patch(isocaps(X,Y,Z,VB2,cutB2),'FaceColor','green','EdgeColor', ...
          'none');
alpha(pB2,0.5);
alpha(qB2,0.5);

% B3
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pB3 = patch(isosurface(X,Y,Z,VB3,cutB3),'FaceColor','green','EdgeColor','none');
qB3 = patch(isocaps(X,Y,Z,VB3,cutB3),'FaceColor','green','EdgeColor', ...
          'none');
alpha(pB3,0.25);
alpha(qB3,0.25);

% B4
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pB4 = patch(isosurface(X,Y,Z,VB4,cutB4),'FaceColor','yellow','EdgeColor','none');
qB4 = patch(isocaps(X,Y,Z,VB4,cutB4),'FaceColor','yellow','EdgeColor', ...
          'none');
alpha(pB4,0.25);
alpha(qB4,0.25);




% A    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pA = patch(isosurface(X,Y,Z,VA,cutA),'FaceColor','red','EdgeColor','none');
qA = patch(isocaps(X,Y,Z,VA,cutA),'FaceColor','red','EdgeColor', ...
         'none');
alpha(pA,0.5);
alpha(qA,0.5);

% C
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pC = patch(isosurface(X,Y,Z,VC,cutC),'FaceColor','blue','EdgeColor','none');
qC = patch(isocaps(X,Y,Z,VC,cutC),'FaceColor','blue','EdgeColor', ...
          'none');
alpha(pC,0.5);
alpha(qC,0.5);


