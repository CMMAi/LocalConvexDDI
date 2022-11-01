function [mech_stress] = updateMechanicalStress(obj, optim_stress, eta)

% Update mechanical stresses.

arguments
    obj
    optim_stress (3,:,:,:) double
    eta (:,:) double
end

B = obj.BMatrices;
C = obj.CMatrix;
weight = obj.weights;


[~, ngp, numberElements,numLoadSteps] = size(optim_stress);

mech_stress = zeros(3,ngp,numberElements,numLoadSteps);

for h = 1:numLoadSteps
    for e = 1:numberElements
        for q = 1:ngp

            stressIn = optim_stress(:,q,e,h);
            w = weight(q,e);
            BMat = B(:,:,q,e);
            eta_temp = eta(:,h);

            [stressOut] = adjust_stress(stressIn, w, C, BMat, eta_temp);
            mech_stress(:,q,e,h) = stressOut;
            
        end
    end
end

function [stressOut] = adjust_stress(stressIn, weight, CTan, BMat, eta)
% Adjust value of stress

LHS = weight * (CTan^-1);
RHS = weight * (CTan^-1)*stressIn + weight * BMat * eta;
stressOut = LHS \ RHS;

end

end
