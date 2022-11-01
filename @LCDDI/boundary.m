function boundary(obj, force, prescribedDof)

% Read boundary conditions
arguments
    obj
    force (:,:) double % (nodal forces, load steps)
    prescribedDof double {mustBeVector(prescribedDof)}
end

obj.force = force;
obj.prescribedDof = prescribedDof;

obj.numLoadSteps = size(obj.force, 2);

end