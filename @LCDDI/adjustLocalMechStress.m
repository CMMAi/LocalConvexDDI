function [mech_stress] = adjustLocalMechStress(obj, optim_stress, eta)

% Adjust mechanical stresses from optimal stresses.

arguments
    obj
    optim_stress (3,:,:,:) double
    eta (:,:) double
end

ngp = obj.ngp;
numberElements = obj.numberElements;
numLoadSteps = obj.numLoadSteps;

B = obj.BMatrices;
C = obj.CMatrix;
weight = obj.weights;


mech_stress = zeros(size(optim_stress));

for h = 1:numLoadSteps
    for e = 1:numberElements
        for q = 1:ngp
            w = weight(q,e);

            A = w*(C^-1);
            b = ( w*(C^-1) * optim_stress(:,q,e,h) + ...
                w * B(:,:,q,e) * eta(:,h) );
            
            % updated mechanical stress
            mech_stress(:,q,e,h) = A\b;
        end
    end
end

end