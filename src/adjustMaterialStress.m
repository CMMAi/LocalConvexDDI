function [mat_stress] = adjustMaterialStress(numMatStates, mech_stress, ...
    mapping, C, weight)
%  Update material stresses based on given mechanical stresses.

arguments
    numMatStates (1,1) double {mustBeInteger(numMatStates)}
    mech_stress (3,:,:) double
    mapping (:,:) double
    C (3,3) double
    weight double {mustBeVector(weight)}
end

mat_stress = zeros(3,numMatStates);



for i = 1:numMatStates
    [elem, loadcase] = find( mapping == i );
    
    A = zeros(3,3); b = zeros(3,1);
    for jj = 1:length(elem)

            A = A + weight(elem(jj))*(C^-1);
            b = b + weight(elem(jj))*(C^-1)*mech_stress(:,elem(jj),loadcase(jj));
    end

    % updated material stress
    mat_stress(:,i) = A\b; 
end

end