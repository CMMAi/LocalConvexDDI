function [mat_stress, eta, optim_strain, optim_stress, LCDD_weight,...
    opt_cell, mech_cell, mat_cell, weight_cell, tNbrSearch] = ...
materialStressSolver(GDof, numLoadSteps, C, pseudoK,numberElements, ...
prescribedDof, weight, B, mat_stress, numMatStates, mappingK, ...
force, tol ,NumK, XiBar, MuBar,mech_strain, mech_stress, mat_strain)

% This function solves for Lagrange multipliers, optimal stress & mat stress
% for the LCDDI algorithm

arguments
    GDof (1,1) double {mustBeInteger(GDof)}
    numLoadSteps (1,1) double {mustBeInteger(numLoadSteps)}
    C (3,3) double
    pseudoK (:,:) double
    numberElements (1,1) double {mustBeInteger(numberElements)}
    prescribedDof double {mustBeVector(prescribedDof)}
    weight double {mustBeVector(weight)}
    B (3,:,:) double
    mat_stress (3,:) double
    numMatStates (1,1) double {mustBeInteger(numMatStates)}
    mappingK (:,:,:) double
    force (:,:) double
    tol (1,1) double
    NumK (1,1) double {mustBeInteger(NumK)}
    XiBar (1,1) double
    MuBar (1,1) double
    mech_strain (3,:,:) double
    mech_stress (3,:,:) double
    mat_strain (3,:) double
end


eta = zeros(GDof,numLoadSteps);
old_mat_stress = zeros( size(mat_stress) );


opt_cell = cell(500, 2);
mech_cell = cell(500, 2);
mat_cell = cell(500, 2);
weight_cell_1 = cell(5, 1);
weight_cell_2 = cell(5, 1);



activeDof = setdiff((1:GDof)',prescribedDof);
LHS = pseudoK(activeDof,activeDof); % left-hand side
L = chol(LHS,'lower');

% time for nearest-neighbor-search
tNbrSearch = 0;


delta = 1; k = 0;
while delta > tol
    k = k + 1;
    
    % Apply LCDD
    startNbrSearch =  cputime;

    [ ~, ~, optim_strain, optim_stress,LCDD_weight] =  ...
        localConvDataSearch(NumK, XiBar, MuBar, numberElements, ...
        mech_strain, mech_stress, numLoadSteps, mat_strain, mat_stress, C, mappingK);

    tNbrSearch = tNbrSearch + (cputime - startNbrSearch);


    [term] = forceFromOptimalStress(weight, B, GDof, numberElements, ...
        numLoadSteps, optim_stress);

    % Solve for Lagrange Multipliers
    parfor i = 1:numLoadSteps
        b = force(:,i) - term{i}(:);
        b = b(activeDof);
        
        x = L'\(L\b);
        eta(activeDof,i) = x;
    end
    

    % Adjust mechanical stress from optimal stress
    [mech_stress] = adjustLocalMechStress(optim_stress, eta, ...
        B, C, weight);
    

    % Adjust material stress from mechanical stress
    mapping = mappingK(:, :, 1);
    [mat_stress] = adjustMaterialStress(numMatStates, mech_stress, ...
        mapping, C, weight);
   
    
    delta = max(abs(mat_stress - old_mat_stress), [], 2);
    old_mat_stress = mat_stress;
    if k <= length(weight_cell_1)
        weight_cell_1{k} = LCDD_weight;
    end
    
    % The mechanical & material stress values are considered to be stable 
    % when "delta" is lower than the imposed tolerance
    for c = 1:-1+length(opt_cell)
        opt_cell{c, 1} = opt_cell{c+1, 1};
        mech_cell{c, 1} = mech_cell{c+1, 1};
        mat_cell{c, 1} = mat_cell{c+1, 1};

        opt_cell{c, 2} = opt_cell{c+1, 2};
        mech_cell{c, 2} = mech_cell{c+1, 2};
        mat_cell{c, 2} = mat_cell{c+1, 2};
    end
    for c = 1:-1+length(weight_cell_2)
        weight_cell_2{c} = weight_cell_2{c+1};
    end
    opt_cell{end, 1} = optim_strain;
    mech_cell{end, 1} = mech_strain;
    mat_cell{end, 1} = mat_strain;
    weight_cell_2{end} = LCDD_weight;

    opt_cell{end, 2} = optim_stress;
    mech_cell{end, 2} = mech_stress;
    mat_cell{end, 2} = mat_stress;
end

% delete rows with empty cells
opt_cell = opt_cell(~any(cellfun('isempty', opt_cell), 2), :);
mech_cell = mech_cell(~any(cellfun('isempty', mech_cell), 2), :);
mat_cell = mat_cell(~any(cellfun('isempty', mat_cell), 2), :);

weight_cell = [weight_cell_1; weight_cell_2];

fprintf('---> (mat stress values stable after %d iterations)\n',k);

end
