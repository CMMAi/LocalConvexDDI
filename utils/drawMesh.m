function drawMesh(nodeCoordinates, elementNodes, type, format)
% Draw mesh discretized with 'T3', 'T6' or 'Q4' elements.

arguments
    nodeCoordinates (:,2) double
    elementNodes (:,:) double
    type char
    format char = 'k-'
end


switch type
    case {'T3','T6'}
        seg1 = [1,2,3,1];
    case 'Q4'
        seg1 = [1,2,3,4,1];
    otherwise
		disp('Type is not supported yet')
end

for e = 1:length(elementNodes(:,1))

    X = nodeCoordinates(elementNodes(e,seg1), 1);
    Y = nodeCoordinates(elementNodes(e,seg1), 2);

    plot(X, Y, format)

    hold on
end

axis equal