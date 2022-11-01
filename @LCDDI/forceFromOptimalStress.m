function [term] = forceFromOptimalStress(obj, optim_stress)

% This function computes the term moved to the right-hand side 
% in modified version of Eq25 (DDI)

arguments
    obj
    optim_stress (3,:,:,:) double
end

numberElements = obj.numberElements;
GDof = obj.GDof;
numLoadSteps = obj.numLoadSteps;

weight = obj.weights;
B = obj.BMatrices;

ngp = obj.ngp;



term = cell(numLoadSteps,1); % term moved to right-hand side in modified version of Eq25 (DDI)

for h = 1:numLoadSteps

    local_term = zeros(GDof,1);
    for e = 1:numberElements

        for q = 1:ngp
        
            local_term(:) = local_term(:) + ...
                ( weight(q,e)*B(:,:,q,e) )' * optim_stress(:,q,e,h);
        end
    end
    term{h} = local_term;
end

end % end function
