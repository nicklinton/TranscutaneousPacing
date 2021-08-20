clear,close all;clc

isDebug = true;

% sort out the torso as a triangulation object.
load('ThoraxAndHeartTriangulations.mat');
load('ThoraxRegions.mat')

% for debugging
if false
    triHeart = trisphere(3);
end
% heartmodel
    trisurf(triHeart);
    axis equal
    hold on
    hTorso = trisurf(triTorso,(ElementID2RegionID-1)*127+1,'FaceAlpha',0.5)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create torso pde model

torso_model = createpde();
% Geometry
    geometryFromMesh(torso_model,triTorso.Points',triTorso.ConnectivityList',ElementID2RegionID);
    generateMesh(torso_model); %clean up the mesh
    

    %tetramesh(triTorsoMesh, 'FaceColor','none')
    triTorsoMesh = triangulation(torso_model.Mesh.Elements(1:4,:)',torso_model.Mesh.Nodes');
    skinFaces = freeBoundary(triTorsoMesh);
    P = triTorsoMesh.Points;
    figure
    trisurf(skinFaces,P(:,1),P(:,2),P(:,3),'FaceAlpha',0.5);
    hold on; axis equal
    trisurf(triHeart);
    


   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create pde model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% PDESystemSize
    % = 1 (which is the default)

% Geometry

   
% first identify a group (='cell') of elements whose surface forms a region

msh = torso_model.Mesh;
p = [132.645 , -320.86 , -196.4];
p_1 = [-176.208,-164.018,-211.69];
% p_2 = [....]
distMax = 30;
distMax_1 = 80;

nodes = msh.Nodes;
elements = msh.Elements;
elemCoords = cell(1,3);
elemCoordsGeometricCenter = zeros(size(elements,2),3);
for i = 1:3 % quivalent to i = x,y,z
    % Find the coordinates of the geometric centers of all elements of the
    % mesh.
    elemCoords{i} = reshape(nodes(i,elements),10,[]);
    % Compute the mean of each column of this array to get a vector of the
    % x-coordinates of the element geometric centers.
    elemCoordsGeometricCenter(:,i) = mean(elemCoords{i});
end

% Assign region ID 1 to all ...
ElementIdToRegionId = ones(1,size(elements,2));


vec = bsxfun(@minus,p, elemCoordsGeometricCenter);
vec_1 = bsxfun(@minus,p_1, elemCoordsGeometricCenter);
distSq = sum(vec.*vec,2);
distSq_1= sum(vec_1.*vec_1,2);
isNear = distSq < (distMax*distMax);
isNear_1 = distSq_1 < (distMax_1*distMax_1);
ElementIdToRegionId(isNear) = 2;
ElementIdToRegionId(isNear_1) = 2;
torso_model2 = createpde();

t = geometryFromMesh(torso_model2,nodes,elements,ElementIdToRegionId);
%generateMesh(torso_model); %clean up the mesh
%figure
% pdegplot(torso_model2,'CellLabels','on','FaceLabels','on','FaceAlpha',0.5)
% FaceID = cellFaces(torso_model2.Geometry,2,'external');
heart_torso_model = createpde();
heart_torso_geometry = addCell(t,h);
heart_torso_model.Geometry = heart_torso_geometry;
pdegplot(heart_torso_model,'CellLabels','on','FaceAlpha',0.4)%'FaceLabels','on','FaceAlpha',0.4);
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply boundary condition
% EquationCoeeficients-LaplaceEquation

specifyCoefficients(heart_torso_model,'m',0,'d',0,'c',1,'a',0,'f',0,'Cell',2);
specifyCoefficients(heart_torso_model,'m',0,'d',0,'c',1,'a',0,'f',0,'Cell',1);
specifyCoefficients(heart_torso_model,'m',0,'d',0,'c',0.1,'a',0,'f',0,'Cell',3);
% IsTimeDependent - 0 *this is the default)

% BoundaryConditions
%bottom
applyBoundaryCondition(heart_torso_model, ...
    'dirichlet', 'Face', 99, 'u', 1000);
% top:
applyBoundaryCondition(heart_torso_model, ...
    'dirichlet', 'Face', 98, 'u', 0);
% left & right:
applyBoundaryCondition(heart_torso_model, ...
    'neumann', 'Face', 1, 'q', 0.000, 'g', 0.000);
% InitialConditions ? just put zero everywhere other than boundary
% put in code
% Mesh for solution-generate mesh
mesh = generateMesh(heart_torso_model);%torso_model2,'CellLabels','on','FaceLabels','on','FaceAlpha',0.5)

% solve pde
results = solvepde(heart_torso_model);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cut in slice (gradient plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = meshgrid(-300:25:200,-500:25:200,-650:25:100);
%save xyz.mat = 
V = interpolateSolution(results,X,Y,Z);
V = reshape(V,size(Z)); 
%countour slice

colormap jet
contourslice(X,Y,Z,V,[],[],-250:15:-130)
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
[gradx,grady,gradz] = evaluateCGradient(results,X,Y,Z);


gradx = reshape(gradx,size(X));
grady = reshape(grady,size(Y));
gradz = reshape(gradz,size(Z));

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