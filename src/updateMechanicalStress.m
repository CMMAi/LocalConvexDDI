function [mech_stress] = updateMechanicalStress(...
    optim_stress, weight, B, eta, C, numberElements, numLoadSteps)
% Compute mechanical stresses.


arguments
    optim_stress (3,:,:) double
    weight double {mustBeVector(weight)}
    B (3,:,:) double
    eta (:,:) double
    C (3,3) double
    numberElements (1,1) double {mustBeInteger(numberElements)}
    numLoadSteps (1,1) double {mustBeInteger(numLoadSteps)}
end



mech_stress = zeros(3,numberElements,numLoadSteps);

for h = 1:numLoadSteps
    for e = 1:numberElements
        mech_stress(:,e,h) = weight(e)*(C^-1)\( weight(e)*(C^-1)*...
            optim_stress(:,e,h) + weight(e)*B(:,:,e)*eta(:,h) );
    end
end

end % end function