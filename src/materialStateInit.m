function [mat_strain, mat_stress, mappingID] = materialStateInit(mech_strain, ...
    numberElements, numLoadSteps, numMatStates)

% This function initializes the material states (strain & stress) as well
% as the mapping between them and the mechanical states.

arguments
    mech_strain (3,:,:) double
    numberElements (1,1) double {mustBeInteger(numberElements)}
    numLoadSteps (1,1) double {mustBeInteger(numLoadSteps)}
    numMatStates (1,1) double {mustBeInteger(numMatStates)}
end


% putting all strain data in a vector in order to apply kmeans
all_strains = zeros(numberElements*numLoadSteps, 3);

k = 0;
for e = 1:numberElements
    for j = 1:numLoadSteps
        k = k+1;
        all_strains(k,:)= mech_strain(:,e,j)';
    end
end

% material strains & initial clustering
[cluster, center]= kmeans(all_strains,numMatStates,'MaxIter',1000);

mat_strain = center';

cluster = reshape(cluster,numLoadSteps,[]); % reshaping vector into a 2D array
mappingID = cluster';

% stresses initialized to 0s
mat_stress = zeros(size(mat_strain));

end % end function