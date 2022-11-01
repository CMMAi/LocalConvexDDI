
clear; clc; close all;

UTILS_PATH = './utils/';
addpath(UTILS_PATH)


% location for saving results
RESULTS_FOLDER = './results/';

rng(1);

%% INPUT DATA

MAT_FILE = './data/dogBone-elastoplastic.mat';
data = load(MAT_FILE);
elemType = 'T3';

nodeCoordinates = data.nodeCoordinates;
elementNodes = data.elementNodes;
thickness = data.thickness;

prescribedDof = data.prescribedDof;
force = data.force;

mech_strain = data.mech_strain;


figure(1)
drawMesh(data.nodeCoordinates, data.elementNodes, elemType)


%%
C = 1e7*eye(3,3);
numMatStates = 100;
NumK = 36;
tol = 5e-2; % convergence criterion for iterative solver

% setup
solver = LCDDI(elemType);
solver.geometry(nodeCoordinates,elementNodes,thickness)
solver.boundary(force, prescribedDof)
solver.strainField(mech_strain);


% run LCDDI
solver.setParams(C, numMatStates, NumK, tol)
[results] = solver.run();



%% SAVE RESULTS
if not(isfolder(RESULTS_FOLDER))
    mkdir( RESULTS_FOLDER );
end

fileName = "database_" + num2str(numMatStates) + "_k" + num2str(NumK);
filePath = fullfile(RESULTS_FOLDER, fileName);

save(filePath, '-v7.3')

