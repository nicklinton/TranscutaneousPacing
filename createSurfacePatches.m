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
    hTorso = trisurf(triTorso,(ElementID2RegionID-1)*127+1,'FaceAlpha',0.5);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create torso pde model

torso_model = createpde();
% Geometry
    geometryFromMesh(torso_model,triTorso.Points',triTorso.ConnectivityList');
    generateMesh(torso_model); %clean up the mesh
    

    %tetramesh(triTorsoMesh, 'FaceColor','none')
    triTorsoMesh = triangulation(torso_model.Mesh.Elements(1:4,:)',torso_model.Mesh.Nodes');
    skinFaces = freeBoundary(triTorsoMesh);
    P = triTorsoMesh.Points;

% now use 
    % userdata = selectregionsonshell(hTorso, 'new',{'frontPatch','backPatch'})
    % selectregionsonshell(hTorso, 'off')
    % selectregionsonshell(hTorso, 'on')
    % selectregionsonshell(hTorso, 'nextregion')
    % userdata = selectregionsonshell(hTorso, 'getdata')

% then convert userdata into the regions
    % ElementID2RegionID = zeros(1,size(skinFaces,1));
    % ElementID2RegionID(userdata.faceRegions(:,1)) = 1;
    % ElementID2RegionID(userdata.faceRegions(:,2)) = 2;
    % ElementID2RegionID = 1+ElementID2RegionID;

% or load SkinFaceID2RegionID
load('SkinFaceID2RegionID_1.mat')

% for each face, find the attached tetrahedron in the original mesh



    figure
    hTorso = trisurf(skinFaces,P(:,1),P(:,2),P(:,3),(SkinFaceID2RegionID-1)*127+1,'FaceAlpha',0.5);
        hold on; axis equal
    trisurf(triHeart);

% now find the mesh elements with a patch surface

frontPatch = skinFaces(SkinFaceID2RegionID==2,:);
frontPatchID = triFaceAttachments(triTorsoMesh, frontPatch);
frontPatchID = cell2mat(frontPatchID);

frontPatchMesh = triangulation(triTorsoMesh.ConnectivityList(frontPatchID,:),triTorsoMesh.Points);

backPatch = skinFaces(SkinFaceID2RegionID==3,:);
backPatchID = triFaceAttachments(triTorsoMesh, backPatch);
backPatchID = cell2mat(backPatchID);

backPatchMesh = triangulation(triTorsoMesh.ConnectivityList(backPatchID,:),triTorsoMesh.Points);

figure
    hTorso = trisurf(skinFaces,P(:,1),P(:,2),P(:,3),(SkinFaceID2RegionID-1)*127+1,'FaceAlpha',0.5);
        hold on; axis equal
        tetramesh(frontPatchMesh, 'FaceColor','b')
        tetramesh(backPatchMesh, 'FaceColor','r')

ElementID2RegionID = ones(size(triTorsoMesh,1),1);
ElementID2RegionID(frontPatchID) = 2;
ElementID2RegionID(backPatchID) = 3;



















