clear all; close all;
load('torsoTetraMesh.mat')
load('ElementID2RegionID_1.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create torso pde model

torso_model = createpde();

% Geometry
    geometryFromMesh(torso_model,triTorsoMesh.Points',triTorsoMesh.ConnectivityList',ElementID2RegionID);
    generateMesh(torso_model); %clean up the mesh
    figure
    pdegplot(torso_model,'FaceLabels','on')

% PDESystemSize
    % = 1 (which is the default)

% EquationCoeeficients-LaplaceEquation
    specifyCoefficients(torso_model,    'm',0,...
                                        'd',0,...
                                        'c',1,...
                                        'a',0,...
                                        'f',0);
% Boundary Conditions
    % Front patch faces ...
    frontPatch = cellFaces(torso_model.Geometry,2,'external');
    applyBoundaryCondition(torso_model,'dirichlet', 'Face', frontPatch, 'u', 100);
    % Back patch faces ...
    backPatch = cellFaces(torso_model.Geometry,3,'external');
    applyBoundaryCondition(torso_model,'dirichlet', 'Face', backPatch, 'u', 0);
    
% solve pde
    results = solvepde(torso_model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

