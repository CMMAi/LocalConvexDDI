function [mat_strain] = updateMaterialStrain(obj, mech_strain, mappingK)

% Update the material strains.

arguments
    obj
    mech_strain (3,:,:,:) double
    mappingK (:,:,:,:) double
end


numMatStates = obj.numMatStates;
weight = obj.weights;
C = obj.CMatrix;


mat_strain = zeros(3,numMatStates);

for i = 1:numMatStates
    
    mapping = mappingK(:,:,:,1);
    indices = find( mapping == i );
    sz = size(mapping);
    [gpts, elem, loadcase] = ind2sub(sz, indices);

    A = zeros(3,3);
    b = zeros(3,1);

    for jj = 1:length(indices)
        gp = gpts(jj);
        e = elem(jj);
        h = loadcase(jj);

        w = weight(gp,e);

        A = A + w * C;
        b = b + w * C * mech_strain(:,gp,e,h);
    end
    
    % updated material strain
    mat_strain(:,i) = A\b; 
end

end