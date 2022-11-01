function [mat_stress, eta, optim_strain, optim_stress, LCDD_weight,...
    solvIter] = materialStressSolver(obj, mat_strain, mat_stress, ...
    mech_strain, mech_stress, ...
    mappingK)

% This function solves for material stresses, Lagrange multipliers, 
% as well as optimal stresses

arguments
    obj
    mat_strain (3,:) double
    mat_stress (3,:) double
    mech_strain (3,:,:,:) double
    mech_stress (3,:,:,:) double
    mappingK (:,:,:,:) double
end

GDof = obj.GDof;
numLoadSteps = obj.numLoadSteps;

prescribedDof = obj.prescribedDof;
force = obj.force;

pseudoK = obj.pseudoK;

tol = obj.tol;




eta = zeros(GDof, numLoadSteps);
old_mat_stress = zeros( size(mat_stress) );


activeDof = setdiff((1:GDof)',prescribedDof);
% left-hand side
LHS = pseudoK(activeDof,activeDof);
L = chol(LHS,'lower');


delta = 1e+6;
solvIter = 0;

while delta > tol
    solvIter = solvIter + 1;
    
    % Apply LCDD
    [ ~, ~, optim_strain, optim_stress, LCDD_weight] =  ...
        obj.localConvDataSearch(mech_strain, mech_stress, ...
        mat_strain, mat_stress, mappingK);


    [term] = obj.forceFromOptimalStress(optim_stress);

    % Solve for Lagrange Multipliers
    parfor i = 1:numLoadSteps
        b = force(:,i) - term{i}(:);
        b = b(activeDof);
        
        x = L'\(L\b);
        eta(activeDof,i) = x;
    end
    

    % Adjust mechanical stress from optimal stress
    [mech_stress] = obj.adjustLocalMechStress(optim_stress, eta);
    

    % Adjust material stress from mechanical stress
    mapping = mappingK(:, :, :, 1);
    [mat_stress] = obj.adjustMaterialStress(mech_stress, mapping);
   
    
    delta = max(abs(mat_stress - old_mat_stress), [], 2);
    old_mat_stress = mat_stress;

end

fprintf('---> (mat stress values stable after %d iterations)\n',solvIter);

end
