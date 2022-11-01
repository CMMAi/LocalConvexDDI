function [mat_stress] = adjustMaterialStress(obj, mech_stress, mapping)

%  Update material stresses based on given mechanical stresses.

arguments
    obj
    mech_stress (3,:,:,:) double
    mapping (:,:,:) double
end

numMatStates = obj.numMatStates;
weight = obj.weights;
C = obj.CMatrix;



mat_stress = zeros(3,numMatStates);

for i = 1:numMatStates

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

        A = A + w * (C^-1);
        b = b + w * (C^-1) * mech_stress(:,gp,e,h);
    end

    % updated material stress
    mat_stress(:,i) = A\b; 
end

end