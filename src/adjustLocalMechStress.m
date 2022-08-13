function [mech_stress] = adjustLocalMechStress(optim_stress, eta, ...
    B, C, weight)

% Adjust mechanical stresses from optimal stresses.

arguments
    optim_stress (:,:,:) double
    eta (:,:) double
    B (3,:,:) double
    C (3,3) double
    weight double {mustBeVector(weight)}
end

mech_stress = zeros(size(optim_stress));
[~, numberElements, numLoadSteps]  = size(mech_stress);

for h = 1:numLoadSteps
    for e = 1:numberElements
        
        A = weight(e)*(C^-1);
        b = ( weight(e)*(C^-1) * optim_stress(:,e,h) + ...
            weight(e) * B(:,:,e) * eta(:,h) );
        
        % updated mechanical stress
        mech_stress(:,e,h) = A\b;
    end
end

end