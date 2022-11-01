function strainField(obj, mech_strain)

% Read strain field measurements

arguments
    obj
    mech_strain (:,:,:,:) double % (3, # of gauss pts, # of elements, # of load steps)
end

obj.mech_strain = mech_strain;

end
