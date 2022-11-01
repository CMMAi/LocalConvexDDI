function geometry(obj, nodeCoordinates, elementNodes, thickness)

% Read problem geometry

arguments
    obj
    nodeCoordinates (:,2) double % (nodes, xy)
    elementNodes (:,:) % (elements, nodes)
    thickness (1,1) double
end

msg = "Check the compatibility of 'elementNodes' with the chosen element type.";
if strcmp(obj.elemType, 'T3')
    assert( size(elementNodes,2) == 3, msg );
elseif strcmp(obj.elemType, 'T6')
    assert( size(elementNodes,2) == 6, msg );
end

obj.nodeCoordinates = nodeCoordinates;
obj.elementNodes = elementNodes;
obj.thickness = thickness;

obj.numberNodes = size(obj.nodeCoordinates, 1);
obj.GDof = obj.numberNodes * 2;
obj.numberElements = size(obj.elementNodes, 1);

end

