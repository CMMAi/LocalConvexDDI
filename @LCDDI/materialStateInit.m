function [mat_strain, mat_stress, mappingID] = materialStateInit(obj)

% This function initializes the material states (strain & stress) as well
% as the mapping between them and the mechanical states.


numberElements = obj.numberElements;
numLoadSteps = obj.numLoadSteps;
numMatStates = obj.numMatStates;

ngp = obj.ngp;
mech_strain = obj.mech_strain;



% putting all strain data in a vector in order to apply kmeans
all_strains = zeros(ngp*numberElements*numLoadSteps, 3);

k = 0;
for q = 1:ngp
    for e = 1:numberElements
        for j = 1:numLoadSteps
            k = k+1;
            all_strains(k,:) = mech_strain(:,q,e,j)';
        end
    end
end


% material strains & initial clustering
[cluster, centroid]= kmeans(all_strains,numMatStates,'MaxIter',1000);

% reshaping
mappingID = reshape(cluster, [ngp,numberElements,numLoadSteps]);

mat_strain = centroid';

% stresses initialized to 0s
mat_stress = zeros(size(mat_strain));

end % end function
