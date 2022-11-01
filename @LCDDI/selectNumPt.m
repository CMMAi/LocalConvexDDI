function [ngp] = selectNumPt(elemType)

switch elemType
    case 'T3'
        ngp = 1;
    case 'T6'
        ngp = 3;
end

end

