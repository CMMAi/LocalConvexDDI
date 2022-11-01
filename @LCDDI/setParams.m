function setParams(obj, CMatrix, numMatStates, NumK, tol, options)

% Set parameters for building material database

arguments
    obj
    CMatrix (3,3) double
    numMatStates (1,1) double {mustBeInteger(numMatStates)}
    NumK (1,1) double {mustBeInteger(NumK)}
    tol (1,1) double

    options.XiBar (1,1) double = 1e4;
    options.MuBar (1,1) double = 1e-4;
end

% Check
msg = "Check input parameters:" + newline + ...
    "The # of nearest neighbors must be less than the # of material states.";
assert(NumK < numMatStates, msg)


obj.CMatrix = CMatrix;
obj.numMatStates = numMatStates;
obj.NumK = NumK;
obj.tol = tol;

obj.optConst = struct( ...
    'XiBar', options.XiBar, ...
    'MuBar', options.MuBar ...
    );

end
