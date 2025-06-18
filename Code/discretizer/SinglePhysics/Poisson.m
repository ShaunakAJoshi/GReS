classdef Poisson < SinglePhysics
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
    analSol        
  end

  properties (Constant)
    field = 'Poisson'
  end

  methods
    function obj = Poisson(symmod,params,dofManager,grid,mat,state)
      obj@SinglePhysics(symmod,params,dofManager,grid,mat,state)
    end

    function computeMat(obj)
      % classical
        % general sparse assembly loop over elements for Poromechanics
      subCells = obj.dofm.getFieldCells(obj.field);
      n = sum(obj.mesh.cellNumVerts(subCells).^2);
      Ndof = obj.dofm.getNumDoF(obj.field);
      asbJ = assembler(n,@(el) computeLocalMatrix(obj,el),Ndof,Ndof);
      % loop over cells
      for el = subCells'
        % get dof id and local matrix
        asbJ.localAssembly(el);
      end
      % populate stiffness matrix
      obj.J = asbJ.sparseAssembly();
    end

    function [dofr,dofc,matLoc] = computeLocalMatrix(obj,elID)
      vtkId = obj.mesh.cellVTKType(elID);
      elem = getElement(obj.elements,vtkId);
      [gradN,dJW] = getDerBasisFAndDet(elem,elID,1);
      matLoc = pagemtimes(gradN,'ctranspose',gradN,'none');
      matLoc = matLoc.*reshape(dJW,1,1,[]);
      matLoc = sum(matLoc,3);
      % get global DoF
      nodes = obj.mesh.cells(elID,1:obj.mesh.cellNumVerts(elID));
      dof = obj.dofm.getLocalDoF(nodes,obj.fldId);
      dofr = dof; dofc = dof;
    end

    function computeRhs(obj)
      obj.rhs = obj.J*obj.state.data.u(ents);
    end

    function setState(obj)
      % add poromechanics fields to state structure
      obj.state.data.u = zeros(obj.mesh.nNodes,1);
      obj.state.data.absErr = zeros(obj.mesh.nNodes,1);
    end

    function updateState(obj,dSol)
      % Update state structure with last solution increment
      ents = obj.dofm.getActiveEnts(obj.field);
      obj.state.data.var(ents) = obj.state.data.var(ents) + dSol(getDoF(obj.dofm,obj.field));
      if ~isempty(obj.analSol)
      obj.state.data.absErr(ents) = abs(obj.state.data.var(ents) - obj.analSol(ents));
      end
    end

    function [dof,vals] = getBC(obj,bc,id,t,~)
      switch bc.getCond(id)
        case 'NodeBC'
          ents = bc.getEntities(id);
        otherwise
          error('BC type %s is not available for %s field',cond,obj.field);
      end
      % map entities dof to local dof numbering
      dof = obj.dofm.getLocalEnts(ents,obj.fldId);
      dof = bc.getCompEntities(id,dof);
      vals = obj.getBCVals(bc,id,t);
    end

    function [cellData,pointData] = printState(obj,sOld,sNew,t)
      % append state variable to output structure
      switch nargin
        case 2
          var = sOld.data.u;
        case 4
          % linearly interpolate state variables containing print time
          fac = (t - sOld.t)/(sNew.t - sOld.t);
          var = sNew.data.u*fac+sOld.data.u*(1-fac);
        otherwise
          error('Wrong number of input arguments');
      end
      [cellData,pointData] = obj.buildPrintStruct(var);
    end

    function [cellStr,pointStr] = buildPrintStruct(var)
      nPointData = 1;
      if ~isempty(obj.analSol)
        nPointData = nPointData + 1;
      end
      pointStr = repmat(struct('name', 1, 'data', 1), nPointData, 1);
      cellStr = [];
      % Displacement
      pointStr(1).name = 'u';
      pointStr(1).data = var;
      if ~isempty(obj.analSol)
        pointStr(2).name = 'abs_error';
        pointStr(2).data = abs(var-obj.analSol);
      end
    end

    function setAnalSolution(obj,f)
      c = obj.mesh.coordinates;
      obj.analSol = arrayfun(@(i) f(c(i,1),c(i,2),c(i,3)),1:obj.mesh.nNodes);
    end

    function [L2err,H1err] = computeError(obj)
      assert(~isempty(obj.analSol),['Missing analytical solution for ' ...
        'Poisson model \n'])
      L2err = 0;
      H1err = 0;
      for el = 1:obj.mesh.nCells
        vtkId = obj.mesh.cellVTKType(el);
        elem = getElement(obj.elements,vtkId);
        N = getBasisFinGPoints(el);
        [gradN,dJW] = getDerBasisFAndDet(elem,el,1);
        dofId = obj.mesh.cells(el,:);
        locErr = obj.state.data.absErr(dofId);
        N_err = N*locErr;
        N_err = (sum(N_err.*dJW))^2;
        gradN_err = pagemtimes(gradN,'ctranspose',locErr,'none');
        gradN_err = gradN_err.*reshape(dJW,1,1,[]);
        gradN_err = sum(gradN_err,3);
        gradN_err = norm(gradN_err,2);
        L2err = L2err + N_err;
        H1err = N_err + gradN_err;
      end
      L2err = sqrt(L2err);
      H1err = sqrt(H1err);
    end
  end

  methods (Static)
    function out = getField()
      out = Poisson.field;
    end
  end
end

