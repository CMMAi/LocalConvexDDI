function [pseudo_stiffness] = pseudoStiffness(obj)

% This function computes the pseudo-stiffness matrix for triangular elements

elementNodes = obj.elementNodes;
nodeCoordinates = obj.nodeCoordinates;
thickness = obj.thickness;

numberElements = obj.numberElements;
GDof = obj.GDof;

C = obj.CMatrix;

ngp = obj.ngp;
shapefunc = obj.shapefunc;


[gaussWeights, gaussLocations] = obj.gauss2dTri(num2str(ngp));

pseudo_stiffness = zeros(GDof, GDof);
for e = 1:numberElements                           
    numNodePerElement = length(elementNodes(e,:));
    numEDOF = 2*numNodePerElement;
    elementDof = zeros(1,numEDOF);
    
    for i = 1:numNodePerElement
        elementDof(2*i-1) = 2*elementNodes(e,i)-1;
        elementDof(2*i) = 2*elementNodes(e,i);   
    end

    ke = zeros(numEDOF,numEDOF);
    coord = nodeCoordinates(elementNodes(e,:), :);
    for ig = 1:ngp
        wi = gaussWeights(ig);

        gaussPoint = gaussLocations(ig,:);
        xi2 = gaussPoint(2);
        xi3 = gaussPoint(3);
        
        [~, naturalDerivatives] = shapefunc(xi2,xi3);
        [Jacob, ~, XYderivatives] = obj.Jacobian(coord,naturalDerivatives);

        % B matrix
        B = zeros(3,numEDOF);
        B(1,1:2:numEDOF) = XYderivatives(1,:); 
        B(2,2:2:numEDOF) = XYderivatives(2,:);
        B(3,1:2:numEDOF) = XYderivatives(2,:);
        B(3,2:2:numEDOF) = XYderivatives(1,:);

        % Integrate stiffness matrix
        ke = ke + B' * C * B * thickness * wi * det(Jacob); 
    end

    % pseudo-stiffness matrix
    pseudo_stiffness(elementDof,elementDof) = ...
        pseudo_stiffness(elementDof,elementDof) + ke;
 
end

end % end function
