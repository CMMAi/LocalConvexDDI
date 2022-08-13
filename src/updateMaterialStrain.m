function [mat_strain] = updateMaterialStrain(mech_strain, ...
    numMatStates, mappingK, weight, C)
% Update the material strains.


arguments
    mech_strain (3,:,:) double
    numMatStates (1,1) double {mustBeInteger(numMatStates)}
    mappingK (:,:,:) double
    weight double {mustBeVector(weight)}
    C (3,3) double
end



mat_strain = zeros(3,numMatStates);

for i = 1:numMatStates
    [elem, loadcase] = find( mappingK(:,:,1) == i );

    A = zeros(3,3); b = zeros(3,1);
    for jj = 1:length(elem)
            
            A = A + weight(elem(jj))*C;
            b = b + weight(elem(jj))*C*mech_strain(:,elem(jj),loadcase(jj));
    end

    mat_strain(:,i) = A\b; % updated material strain states
end

end