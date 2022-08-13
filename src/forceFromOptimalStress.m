function [term] = forceFromOptimalStress(weight, B, GDof, numberElements, ...
    numLoadSteps, optim_stress)

% This function computes the term moved to the right-hand side 
% in modified version of Eq25 (DDI)

arguments
    weight (:,1) double
    B (3,:,:) double
    GDof (1,1) double
    numberElements (3,:,:) double
    numLoadSteps (1,1) double {mustBeInteger(numLoadSteps)}
    optim_stress (3,:,:) double
end


term = cell(numLoadSteps,1); % term moved to right-hand side in modified version of Eq25 (DDI)

for h=1:numLoadSteps

    local_term = zeros(GDof,1);
    for e = 1:numberElements
        
        local_term(:) = local_term(:) + ...
            ( weight(e)*B(:,:,e) )' * optim_stress(:,e,h);
    end
    term{h} = local_term;
end

end % end function