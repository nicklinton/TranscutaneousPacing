close all
results;

triTorsoMesh = triangulation(results.Mesh.Elements(1:4,:)',results.Mesh.Nodes');

skinFaces = freeBoundary(triTorsoMesh);
p = triTorsoMesh.Points;

hSkin = trimesh(skinFaces, p(:,1), p(:,2), p(:,3), ...
    'EdgeColor','k' , 'EdgeAlpha',0.2, 'FaceColor',[.5 .5 .5], 'FaceAlpha',0.2);
axis vis3d equal
hold on

load('ThoraxAndHeartTriangulations.mat')
p = triHeart.Points;
vHeart = interpolateSolution(results,p');
hHeart = trisurf(triHeart, 'CData',vHeart);

step = 10;
step2 = 5;
[X,Y,Z] = meshgrid(-200:step:200,-350:step:0,-275:step2:-175);
[X,Y,Z] = meshgrid(-200:step:200,-350:step:0,-225+(step2*(-1:0.2:1)));

V = interpolateSolution(results,X,Y,Z);
V = reshape(V,size(Z));
contourslice(X,Y,Z,V,[],[],-225+(step2*(-1:0.2:1)),0:2.5:60)

[X,Y,Z] = meshgrid(-200:step:200,-350:step:0,-225);
[gradx,grady,gradz] = evaluateGradient(results,X,Y,Z);
gradx = reshape(gradx,size(X));
grady = reshape(grady,size(Y));
gradz = reshape(gradz,size(Z));

k = -200;
quiver3(X,Y,Z,k*gradx,k*grady,k*gradz,5,'r')

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cut in slice (gradient plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = 10;
[X,Y,Z] = meshgrid(-300:10:200,0:10:350,-650:10:100);
%save xyz.mat = 
V = interpolateSolution(results,X,Y,Z);
V = reshape(V,size(Z)); 
%countour slice

colormap jet
contourslice(X,Y,Z,V,[],[],-250:15:-130, 0:10:100)
xlabel('x')
ylabel('y')
zlabel('z')
colorbar
view(-11,14);
title('Contour Plot of Estimated Voltage');
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%gradient plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure;
k = -10;
quiver3(X,Y,Z,k*gradx,k*grady,k*gradz)
axis equal
xlabel 'x'
ylabel 'y'
zlabel 'z'
view(-11,14);
title('Quiver Plot of Estimated Current')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cut in slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y] = meshgrid(-300:5:200,-500:5:200);
Z = -177+0*X+0*Y;
V = interpolateSolution(results,X,Y,Z);
V = reshape(V,size(Y));
figure
surf(X,Y,V,'LineStyle','none');
axis equal
view(0,90)
title('Colored Plot on Tilted Plane')
xlabel('x')
ylabel('y')
colorbar
%slice(X,Y,Z,,[],[],-177,'nearest')