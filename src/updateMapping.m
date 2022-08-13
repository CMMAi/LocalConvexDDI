function [mapping1, mappingK, Total_Dist] = ...
    updateMapping(mat_strain, mat_stress,...
    mech_strain, mech_stress, C, numberElements, numLoadSteps, NumK)

arguments
    mat_strain (3,:) double
    mat_stress (3,:) double
    mech_strain (3,:,:) double
    mech_stress (3,:,:) double
    C (3,3) double
    numberElements (1,1) double {mustBeInteger(numberElements)}
    numLoadSteps (1,1) double {mustBeInteger(numLoadSteps)}
    NumK (1,1) double {mustBeInteger(NumK)}
end

mapping1 = zeros(numberElements, numLoadSteps);
mappingK = zeros(numberElements, numLoadSteps, NumK);

% total distance
Total_Dist = 0; 

Mat_eps = C;
Mat_sig = C^-1;
MatDataSet_E = [mat_strain', mat_stress'];

for h = 1:numLoadSteps
    for e = 1:numberElements
        
        X = MatDataSet_E;
        Y = [mech_strain(1,e,h), mech_strain(2,e,h), ...
             mech_strain(3,e,h), mech_stress(1,e,h), ...
             mech_stress(2,e,h), mech_stress(3,e,h)]';
        
        MBar = zeros(6,6);
        MBar(1:3, 1:3) = Mat_eps;
        MBar(4:6, 4:6) = Mat_sig;
        
        [U, S, V] = svd(MBar); 
        MBarSqrRoot = U * sqrt(S) * V';
        NNLS_Vec_z = MBarSqrRoot*Y;
        test = MBarSqrRoot*X';
        
        [Energyidx, distance] = knnsearch(test', NNLS_Vec_z','K', NumK);

        mapping1(e,h) = Energyidx(1); % index of closest material state
        mappingK(e,h,:) = Energyidx;

        Total_Dist = Total_Dist + distance(1);
    end
end


end