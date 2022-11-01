function initialize(obj)

% Initialize necessary variables

arguments
    obj
end

% strain-displacement matrix & element weights
[obj.BMatrices, obj.weights] = obj.strainDispMatrix();

% pseudo stiffness
obj.pseudoK = obj.pseudoStiffness();


% mechanical stress
obj.mech_stress = zeros( size( obj.mech_strain ) );


% initialize material states & mapping
[obj.mat_strain, obj.mat_stress, ~] = obj.materialStateInit();


[obj.mappingK, ~] = obj.updateMapping(obj.mat_strain, obj.mat_stress, ...
    obj.mech_strain, obj.mech_stress);




end

