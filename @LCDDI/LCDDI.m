classdef LCDDI < handle
    % A class for building databases using the "LOCAL-CONVEXITY DATA-DRIVEN
    % IDENTIFICATION" method.

    % PROPERTIES
    properties
        elemType

        nodeCoordinates
        elementNodes
        thickness

        force
        prescribedDof

        mech_strain
    end

    properties (SetAccess=protected)
        numberNodes
        GDof
        numberElements
        numLoadSteps
    end

    properties
        shapefunc
        ngp
    end

    properties
        CMatrix
        numMatStates
        NumK
        tol % tolerance

        optConst % optional constant
    end

    properties (SetAccess=protected)
        pseudoK
        BMatrices
        weights

        mech_stress
        mat_strain
        mat_stress

        mappingK
    end






    % METHODS
    methods
        % Class Constructor
        function obj = LCDDI(elemType)
            
            % Check
            expectedTypes = {'T3','T6'};
            msg = "Only T3 and T6 elements are currently supported";
            assert( any(validatestring(elemType,expectedTypes)), msg )


            obj.elemType = elemType;
            obj.shapefunc = obj.shapeFunction( obj.elemType );
            obj.ngp = obj.selectNumPt( obj.elemType );

        end
    end


    methods
        % Main functions to run algorithm
        geometry(obj, nodeCoordinates, elementNodes, thickness)

        boundary(obj, force, prescribedDof)

        strainField(obj, mech_strain)

        setParams(obj, CMatrix, numMatStates, NumK, tolerance, options)

        [results] = run(obj)

    end


    methods (Access=protected)
        % helper functions for initialization
        initialize(obj)

        [pseudo_stiffness] = pseudoStiffness(obj)

        [mat_strain, mat_stress, mappingID] = materialStateInit(obj)

        [mapping1, mappingK, Total_Dist] = updateMapping(obj, ...
            mat_strain, mat_stress,...
            mech_strain, mech_stress)

        [B, weight] = strainDispMatrix(obj)


    end

    methods (Access=protected)
        [mech_stress] = updateMechanicalStress(obj, optim_stress, eta)

        [mat_strain] = updateMaterialStrain(obj, mech_strain, mappingK)
        
        [mat_stress] = adjustMaterialStress(obj, mech_stress, mapping)


        [mat_stress, eta, optim_strain, optim_stress, LCDD_weight,...
            solvIter] = materialStressSolver(obj, mat_strain, mat_stress, ...
            mech_strain, mech_stress, ...
            mappingK)

        [Local_LCDD_Diff_Max, GlobalEnergyDensity, ...
            M_eps, M_sig, LCDD_weight] =  localConvDataSearch(obj, ...
            mech_strain, mech_stress, ...
            mat_strain, mat_stress, ...
            mappingK)

        [term] = forceFromOptimalStress(obj, optim_stress)

        [mech_stress] = adjustLocalMechStress(obj, optim_stress, eta)

    end

    methods (Static)
        [shapefunc] = shapeFunction( elemType )

        [ngp] = selectNumPt( elemType )
    
        [JacobianMatrix,invJacobian,XYDerivatives] = ...
        Jacobian(nodeCoordinates,naturalDerivatives)

        [weights,locations] = gauss2dTri(option)
    end


%     methods (Static)
%         % function to check convergence
%         [checkFlag, mapping_LCDD, mappingID] = checkConverg(...
%             mappingID, new_mappingID, mapping_LCDD,new_mapping_LCDD)
%     end

end % end class definition

