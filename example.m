
clear; clc; close all;

SRC_PATH = './src/';
addpath(SRC_PATH)

% location for saving results
RESULTS_FOLDER = './results/';

rng(1);

%% INPUT DATA
MAT_FILE = './data/dogBone-elastoplastic.mat';
data = load(MAT_FILE);

geometry = struct( ...
    'nodeCoordinates', data.nodeCoordinates,...
    'numberNodes', data.numberNodes, ...
    'elementNodes', data.elementNodes, ...
    'numberElements', data.numberElements, ...
    'thickness', data.thickness ...
    );

boundaryData = struct( ...
    'prescribedDof', data.prescribedDof, ...
    'force', data.force, ...
    'numLoadSteps', data.numLoadSteps ...
    );

figure(1)
drawMesh(data.nodeCoordinates, data.elementNodes, 'T3')

%% PARAMETERS
numMatStates = 100;
NumK = 1;
C = 1e7*eye(3,3);
tol = 5e-2; % convergence criterion for iterative solver
   
params = struct( ...
    'numMatStates', numMatStates, ...
    'NumK', NumK, ...
    'C', C, ...
    'tolerance', tol ...
    );


% run LCDDI
[results] = runLCDDI(geometry, boundaryData, data.mech_strain, params);


%% SAVE RESULTS
if not(isfolder(RESULTS_FOLDER))
    mkdir( RESULTS_FOLDER );
end

fileName = "database_" + num2str(numMatStates) + "_k" + num2str(NumK);
filePath = fullfile(RESULTS_FOLDER, fileName);

save(filePath, '-v7.3')


