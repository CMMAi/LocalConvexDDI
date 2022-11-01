function [Local_LCDD_Diff_Max, GlobalEnergyDensity, M_eps, M_sig, LCDD_weight] =  ...
    localConvDataSearch(obj, mech_strain, mech_stress, ...
    mat_strain, mat_stress, ...
    mappingK)

% Perform Local Convex Data Search

arguments
    obj
    mech_strain (3,:,:,:) double
    mech_stress (3,:,:,:) double
    mat_strain (3,:) double
    mat_stress (3,:) double
    mappingK (:,:,:,:) double
end

C = obj.CMatrix;

numberElements = obj.numberElements;
numLoadSteps = obj.numLoadSteps;
numMatStates = obj.numMatStates;
NumK = obj.NumK;

ngp = obj.ngp;

XiBar = obj.optConst.('XiBar');
MuBar = obj.optConst.('MuBar');



Mat_eps = C; Mat_sig = C^-1;
M_epsInput = mech_strain; M_sigInput = mech_stress;

M_eps = zeros(size(mech_strain));
M_sig = zeros(size(mech_stress));


LCDD_weight = zeros(ngp,numberElements,numLoadSteps,numMatStates);

GlobalEnergyDensity = 0;
Local_LCDD_Diff_Max = 0;

% transforming material state data into form usable by LCDD algorithm
MatDataSet_E = [mat_strain', mat_stress'];


for h = 1:numLoadSteps % total loading conditions
    for e = 1:numberElements               
        for q = 1:ngp

            X = MatDataSet_E;
            Y = [mech_strain(1,q,e,h), mech_strain(2,q,e,h), ...
                mech_strain(3,q,e,h), mech_stress(1,q,e,h), ...
                mech_stress(2,q,e,h), mech_stress(3,q,e,h)]';
    
            MBar = zeros(6,6);
            MBar(1:3, 1:3) = Mat_eps;
            MBar(4:6, 4:6) = Mat_sig;
    
            [U, S, V] = svd(MBar); 
            MBarSqrRoot = U * sqrt(S) * V';
            NNLS_Vec_z = MBarSqrRoot*Y;
    %         test = MBarSqrRoot*X';
            
            Energyidx = squeeze ( mappingK(q,e,h,:) );
    %         Energyidx = knnsearch(test', NNLS_Vec_z','K', NumK);
        
            knn_S = MatDataSet_E(Energyidx,:)';
            

            if (NumK == 1)
                % No need to solve for NNLS problem
                NNLS_Wts = 1;
    
            else
    
                NNLS_Mtr_Z = MBarSqrRoot*knn_S;
                NNLS_Xi = XiBar.*trace(NNLS_Mtr_Z'*NNLS_Mtr_Z)/NumK;                        
                NNLS_Aug_Vec_z = [NNLS_Vec_z; sqrt(NNLS_Xi)];
                NNLS_Aug_Mtr_Z = [NNLS_Mtr_Z; sqrt(NNLS_Xi).*ones(1, NumK)];
        
        
                options = optimset('Display','off');
                if (MuBar <= 0 ) % NNLS               
        
                    NNLS_Wts = lsqnonneg(NNLS_Aug_Mtr_Z, NNLS_Aug_Vec_z, options);
        
                else % NNLS + Tikhonov regularization  
         
                    NNLS_Mu = MuBar.*trace(NNLS_Mtr_Z'*NNLS_Mtr_Z)./NumK;                
                    NNLS_Aug_Vec_z_TR = [NNLS_Aug_Vec_z; zeros(1, NumK)'];             
                    NNLS_Aug_Mtr_Z_TR = [NNLS_Aug_Mtr_Z; sqrt(NNLS_Mu).*eye(NumK)];
                    NNLS_Wts = lsqnonneg(NNLS_Aug_Mtr_Z_TR, NNLS_Aug_Vec_z_TR);
                end

            end
    
            % Interpolated state
            LCDD_State = knn_S * NNLS_Wts;
            Diff = LCDD_State - [M_epsInput(1,q,e,h), M_epsInput(2,q,e,h), ...
                                 M_epsInput(3,q,e,h), M_sigInput(1,q,e,h), ...
                                 M_sigInput(2,q,e,h), M_sigInput(3,q,e,h)]';
            Diff = norm(MBarSqrRoot*Diff);
    
            if(Diff > Local_LCDD_Diff_Max)
                Local_LCDD_Diff_Max = Diff;
            end    
    
            M_eps(1,q,e,h) = LCDD_State(1); M_sig(1,q,e,h) = LCDD_State(4);
            M_eps(2,q,e,h) = LCDD_State(2); M_sig(2,q,e,h) = LCDD_State(5);
            M_eps(3,q,e,h) = LCDD_State(3); M_sig(3,q,e,h) = LCDD_State(6);
            LCDD_weight(q,e,h,Energyidx) = NNLS_Wts;
        end
    end
end

end
