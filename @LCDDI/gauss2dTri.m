function [weights,locations] = gauss2dTri(option)
% Gauss quadrature in 2D in Triangle element
% option '1'
% option '3' 
% option '4' 
% option '6' 
% locations: Gauss point locations
% weights: Gauss point weights

switch option
    case '6'
    alpha1 = 0.8168475730;
    alpha2 = 0.1081030182;
    beta1 = 0.0915762135;
    beta2 = 0.4459484909;
    gamma3 = 0.1099517437/2;
    gamma4 = 0.2233815897/2;
    locations = [ alpha1 beta1 beta1;
        beta1 alpha1 beta1;
        beta1 beta1 alpha1;
        alpha2 beta2 beta2;
        beta2 alpha2 beta2;
        beta2 beta2 alpha2];
    weights = [ gamma3; gamma3; gamma3; gamma4; gamma4; gamma4;];

    case '4'
    gamma1 = -27/96;
    gamma2 = 25/96;
    locations = [ 1/3 1/3 1/3;
        0.6 0.2 0.2;
        0.2 0.6 0.2;
        0.2 0.2 0.6];
    weights = [gamma1; gamma2; gamma2; gamma2];

    case '3'
    locations = [ 1/2 1/2 0;
        0 1/2 1/2;
        1/2 0 1/2];
    weights = [ 1/6;1/6;1/6]; 

    case '1'
    locations = [1/3 1/3 1/3];
    weights = [1/2];
end

end

