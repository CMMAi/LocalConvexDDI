function [mappingK, Total_Dist] = updateMapping(obj, ...
    mat_strain, mat_stress,...
    mech_strain, mech_stress)

% Update mapping between material states & mechanical states
arguments
    obj
    mat_strain (3,:) double
    mat_stress (3,:) double
    mech_strain (3,:,:,:) double
    mech_stress (3,:,:,:) double
end

numberElements = obj.numberElements;
numLoadSteps = obj.numLoadSteps;

C = obj.CMatrix;
NumK = obj.NumK;

ngp = obj.ngp;



mappingK = zeros( ngp, numberElements, numLoadSteps, NumK );


% total distance
Total_Dist = 0; 

Mat_eps = C;
Mat_sig = C^-1;
MatDataSet_E = [mat_strain', mat_stress'];

for h = 1:numLoadSteps
    for e = 1:numberElements
        for q = 1:ngp
        
            X = MatDataSet_E;
            Y = [mech_strain(1,q,e,h), mech_strain(2,q,e,h), ...
                 mech_strain(3,q,e,h), mech_stress(1,q,e,h), ...
                 mech_stress(2,q,e,h), mech_stress(3,q,e,h)]';
            
            [Energyidx, distance] = local_search(X, Y, Mat_eps, Mat_sig, NumK);
    
            mappingK(q,e,h,:) = Energyidx;
    
            Total_Dist = Total_Dist + distance(1);
        end
    end
end


function [Idx, Dist] = local_search(dataSet, state, Mat_eps, Mat_sig, NumK)
% Find closest point to "state" in "dataSet"

MBar = zeros(6,6);
MBar(1:3, 1:3) = Mat_eps;
MBar(4:6, 4:6) = Mat_sig;

[U, S, V] = svd(MBar); 
MBarSqrRoot = U * sqrt(S) * V';
NNLS_Vec_z = MBarSqrRoot * state;
test = MBarSqrRoot * dataSet';

[Idx, Dist] = knnsearch(test', NNLS_Vec_z','K', NumK);

end


end