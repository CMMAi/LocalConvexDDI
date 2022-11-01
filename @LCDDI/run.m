function [results] = run(obj)

% Run LCDDI algorithm

%% Initialize

mech_strain = obj.mech_strain;

% -------------------------
obj.initialize()

mech_stress = obj.mech_stress;
mat_strain = obj.mat_strain;
mat_stress = obj.mat_stress;

mappingK = obj.mappingK;



%%  Main Body of Solver

results =  struct;
solverIters = [];



iter = 0;
convFlag = false;
while (~convFlag)
    iter = iter + 1;
    fprintf(['Iteration: ',num2str(iter),'\n']);
    


    % Compute material stress, Lagrange multipliers, optimal strains &
    % stress, etc.
    [mat_stress, eta, ~, optim_stress, LCDD_weight, ...
        solvIter] = obj.materialStressSolver(mat_strain, mat_stress,...
        mech_strain, mech_stress, ...
        mappingK);

    solverIters = [solverIters;
                    solvIter];


    % Update mechanical stresses
    [mech_stress] = obj.updateMechanicalStress(optim_stress, eta);


    % Update mapping
    [mappingK_new, ~] = obj.updateMapping(...
        mat_strain, mat_stress, ...
        mech_strain, mech_stress);
    

    % Update material strains
    [mat_strain] = obj.updateMaterialStrain(mech_strain, mappingK_new);


    % Test for convergence (only consider 1 closest point)
    mapA = mappingK(:,:,:,1);
    mapB = mappingK_new(:,:,:,1);

    convFlag = isequal( mapA, mapB );
    if convFlag
        break;
    else
        jj = nnz( mapA ~= mapB );
        fprintf('---> (mapping: %d IDs are different)\n',jj);

        mappingK = mappingK_new;
    end 

end % end Solver


%% Post-Processing

obj.mech_stress = mech_stress;
obj.mat_strain = mat_strain;
obj.mat_stress = mat_stress;

obj.mappingK = mappingK;


% -------------------------
results.('mat_strain') = obj.mat_strain;
results.('mat_stress') = obj.mat_stress;
results.('mech_stress') = obj.mech_stress;

results.('mappingK') = mappingK;
results.('iter') = iter;
results.('solverIters') = solverIters;
results.('LCDD_weight') = LCDD_weight;


end

