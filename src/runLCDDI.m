function [results] = runLCDDI(geometry, boundaryData, mech_strain, params, options)


arguments
    geometry struct 
    boundaryData struct
    mech_strain (3,:,:) double
    params struct

    options.XiBar (1,1) double = 1e4;
    options.MuBar (1,1) double = 1e-4;
end

geomFields = {'nodeCoordinates','numberNodes','elementNodes','numberElements','thickness'};
bcFields = {'prescribedDof','force','numLoadSteps'};
paramFields = {'NumK','C','numMatStates','tolerance'};

inpMsg = 'Please check fieldnames in input structure arrays.';
assert( all(isfield(geometry, geomFields)), inpMsg );
assert( all(isfield(boundaryData, bcFields)), inpMsg );
assert( all(isfield(params, paramFields)), inpMsg );

%% Retrieve Input

% Geometry
nodeCoordinates = geometry.('nodeCoordinates');
numberNodes = geometry.('numberNodes');
GDof = 2 * numberNodes;

elementNodes = geometry.('elementNodes');
numberElements = geometry.('numberElements');
thickness = geometry.('thickness');


% Boundary Conditions
prescribedDof = boundaryData.('prescribedDof');
force = boundaryData.('force');
numLoadSteps = boundaryData.('numLoadSteps');


% Parameters for Solver
NumK = params.('NumK');
C = params.('C');
numMatStates = params.('numMatStates');
tol = params.('tolerance');

XiBar = options.XiBar;
MuBar = options.MuBar;


% Check
msg = "Check input parameters:" + newline + ...
    "The # of nearest neighbors must be less than the # of material states.";
assert(NumK < numMatStates, msg)



%% Initialization

mech_stress = zeros(size(mech_strain));


% Pseudo stiffness
[pseudoK] = pseudoStiffness(...
    GDof,numberElements,elementNodes,nodeCoordinates,C,thickness);

% Initialize material states & mapping
[mat_strain, mat_stress, mapping1] = materialStateInit(...
    mech_strain,numberElements,numLoadSteps,numMatStates);



[~, mappingK, ~] = updateMapping(...
    mat_strain, mat_stress,mech_strain, mech_stress, C, numberElements, numLoadSteps, NumK);


% Strain-displacement matrix & element weights
[B, weight] = strainDispMatrix(GDof, numberElements, elementNodes, ...
    nodeCoordinates, thickness);




%%  Main Body of LCDDI Solver
results =  struct;

startAlgo =  cputime;
timeStressSolver = 0;
timeNbrSearch = 0;



iter = 0;
while true
    iter = iter + 1;
    fprintf(['Iteration: ',num2str(iter),'\n']);
    


    % Compute material stress, Lagrange multipliers, optimal strains &
    % stress, etc.
    
    startStressSolver =  cputime;
    [mat_stress, eta, ~, optim_stress, ~, ...
        ~, ~, ~, ~, tNbrSearch] = ...
        materialStressSolver(...
                            GDof, numLoadSteps, C, pseudoK,numberElements,...
                            prescribedDof, weight, B, mat_stress,...
                            numMatStates, mappingK, force,...
                            tol ,NumK, XiBar, MuBar,mech_strain, ...
                            mech_stress, mat_strain);

    timeStressSolver = timeStressSolver + (cputime - startStressSolver);
    timeNbrSearch = timeNbrSearch + tNbrSearch;

    % Update mechanical stresses
    [mech_stress] = updateMechanicalStress(...
        optim_stress, weight, B, eta, C, numberElements,numLoadSteps);


    % Update mapping
    [mapping1_new, mappingK_new, ~] = updateMapping(...
        mat_strain, mat_stress, mech_strain, mech_stress, C,...
        numberElements, numLoadSteps, NumK);
    

    % Update material strains
    [mat_strain] = updateMaterialStrain(...
        mech_strain, numMatStates, mappingK_new, weight, C);

    % Test for convergence
    if isequal( mappingK(:,:,1), mappingK_new(:,:,1) )
        break;
    else
        jj = nnz( mappingK(:,:,1) ~= mappingK_new(:,:,1) );
        fprintf('---> (mapping: %d IDs are different)\n',jj);

        mapping1 = mapping1_new;
        mappingK = mappingK_new;
    end 

end % end LCDDI Solver

timeAlgo =  cputime - startAlgo;


%% Results for ouput

results.('mat_strain') = mat_strain;
results.('mat_stress') = mat_stress;
results.('mech_stress') = mech_stress;

results.('mapping4LCDD') = mappingK;

results.('timeAlgo') = timeAlgo;
results.('timeStressSolver') = timeStressSolver;
results.('timeNbrSearch') = timeNbrSearch;

results.('iter') = iter;


end


