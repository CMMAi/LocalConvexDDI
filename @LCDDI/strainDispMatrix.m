function [B, weight] = strainDispMatrix(obj)

% This function outputs:
% B: strain-displacement matrix
% weight: element weights


elementNodes = obj.elementNodes;
nodeCoordinates = obj.nodeCoordinates;
GDof = obj.GDof;
numberElements = obj.numberElements;
thickness = obj.thickness;

shapefunc = obj.shapefunc;
ngp = obj.ngp;


[gaussWeights, gaussLocations] = obj.gauss2dTri(num2str(ngp));


B = zeros(3, GDof, ngp, numberElements);
areas = zeros(ngp, numberElements);

for e = 1:numberElements                           
    numNodePerElement = length(elementNodes(e,:));
    numEDOF = 2*numNodePerElement;
    elementDof=zeros(1,numEDOF);

    for i = 1:numNodePerElement
        elementDof(2*i-1)=2*elementNodes(e,i)-1;
        elementDof(2*i)=2*elementNodes(e,i);   
    end
    
    coord = nodeCoordinates(elementNodes(e,:), :);
    for q = 1:ngp
        wi = gaussWeights(q);

        gaussPoint = gaussLocations(q,:);
        xi2 = gaussPoint(2);
        xi3 = gaussPoint(3);
        
        [~, naturalDerivatives] = shapefunc(xi2,xi3);
        [Jacob, ~, XYderivatives] = obj.Jacobian(coord,naturalDerivatives);

        % B matrix
        bee = zeros(3,numEDOF);
        bee(1,1:2:numEDOF) = XYderivatives(1,:); 
        bee(2,2:2:numEDOF) = XYderivatives(2,:);
        bee(3,1:2:numEDOF) = XYderivatives(2,:);
        bee(3,2:2:numEDOF) = XYderivatives(1,:);
    
        % scattering results in array
        B(:,elementDof,q,e) = B(:,elementDof,q,e) + bee;

        areas(q,e) = wi * det(Jacob);
    end
end

weight = areas * thickness;


end