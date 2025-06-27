classdef MeshGlueDual < MeshGlue
  % subclass of Mortar implementing mesh tying between different domains
  % using dual multipliers

  properties (Access = public)
    E
    Jcoupling
  end

  methods

    function obj = MeshGlueDual(id,inputStruct,domains)
      obj@MeshGlue(id,inputStruct,domains);
      assert(isfield(inputStruct,"Physics"), ...
        'Missing Physics field for interface %i',id);
      assert(strcmp(obj.multiplierType,'dual'));
    end
    %
  end

  methods (Access=public)
    %
    function varargout = getJacobian(obj,varargin)
      % get jacobian blocks associated to specific field and specific
      % domain
      % if nargout = 2 -> get master/slave pair of jacobian blocks
      % if nargout = 1 -> return multiplier jacobian
      varargout = cell(1,nargout);
      for i = 1:nargout
        varargout{i} = [];
      end
    end

    function rhs = getRhs(obj,fldId,varargin)
      % return rhs block associated to master/slave field
      switch nargin
        case 2
          % multiplier rhs block
          rhs = [];
        case 3
          domId = varargin{1};
          isDom = ismember(obj.idDomain,domId);
          if all(isDom)
            rhs = obj.rhsMaster{fldId} + obj.rhsSlave{fldId};
          elseif isDom(1)
            rhs = obj.rhsMaster{fldId};
          elseif isDom(2)
            rhs = obj.rhsSlave{fldId};
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end


    function computeMat(obj,~)
      computeMortarMatrices(obj);
      obj.E = computeMortarOperator(obj);
      getCouplingMat(obj);
    end

    function getCouplingMat(obj)
      % return the coupling matrix for static condensation of dual
      % multipliers
      if ~strcmp(obj.multiplierType,'dual')
        return
      end

      solvSlave = obj.solvers(2).getSolver(obj.physics);

      % coupling block between one domain and the other
      obj.Jcoupling = solvSlave.J*obj.E;

      % update master block with condensation term
      solvMaster = obj.solvers(1).getSolver(obj.physics);
      solvMaster.J = solvMaster.J + obj.E'*solvSlave.J*obj.E; 
    end

    function computeRhs(obj)
      % compute rhs contributions for a specified input field
      for i = 1:obj.nFld
        % reset rhs multiplier
        obj.rhsMult{i} = zeros(getNumbMultipliers(obj),1);
        computeRhsMaster(obj,i);
        computeRhsSlave(obj,i);
      end
    end

    function applyBCmaster(obj,bound,bc,t)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(physic),bound,bc,t);
      i = strcmp(obj.physics,physic);
      obj.rhsMaster{i}(bcEnts) = 0;
      obj.Jcoupling(:,bcEnts) = 0;
    end

    function applyBCslave(obj,bound,bc,t)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(physic),bound,bc,t);
      bcEnts = removeSlaveBCdofs(obj,physic,bcEnts);
      i = strcmp(obj.physics,physic);
      obj.rhsSlave{i}(bcEnts) = 0;
      obj.Jcoupling(bcEnts,:) = 0;

      % remove interface slave dofs from matrix system (force zero) 
      dofSlave = obj.mesh.local2glob{2};
      obj.solvers(2).getSolver(obj.physics).applyDirBC([],dofSlave);
    end

    function updateState(obj,du)
      for i = 1:obj.nFld
        % get increment of multipliers from other matrix block
        actMult = getMultiplierDoF(obj);
        obj.multipliers(i).curr(actMult) = obj.multipliers(i).curr(actMult) + x;
      end
    end
  end

  methods (Access = private)

    function initializeJacobianAndRhs(obj)
      [obj.Jmaster,obj.Jslave,obj.Jmult] = deal(cell(obj.nFld,1));
      [obj.rhsMaster,obj.rhsSlave, obj.iniMultipliers] = deal(cell(obj.nFld,1));
      obj.multipliers = repmat(struct('prev',[],'curr',[]),obj.nFld,1);
      obj.totMult = 0;

      for i = 1:obj.nFld
        nDofMaster = getNumDoF(obj.dofm(1),obj.physics(i));
        nDofSlave = getNumDoF(obj.dofm(2),obj.physics(i));
        nDofMult = getNumbMultipliers(obj);
        obj.rhsMaster{i} = zeros(nDofMaster,1);
        obj.rhsSlave{i} = zeros(nDofSlave,1);
        obj.multipliers(i).curr = zeros(nDofMult,1);
        obj.multipliers(i).prev = obj.multipliers(i).curr;
        obj.iniMultipliers{i} = obj.multipliers(i).curr;
        obj.totMult = obj.totMult + nDofMult;
      end
    end

    function computeRhsMaster(obj,i)
      var = getState(obj.solvers(2).getSolver(obj.physics(i)));
      ents = obj.dofm(2).getActiveEnts(obj.physics(i));
      obj.rhsMaster{i} = obj.Jcoupling'*var(ents);
    end

    function computeRhsSlave(obj,i)
      var = getState(obj.solvers(1).getSolver(obj.physics(i)));
      ents = obj.dofm(1).getActiveEnts(obj.physics(i));
      obj.rhsSlave{i} = obj.Jcoupling*var(ents);
      % update rhs of master side with condensation contribution
      obj.rhsMaster{i} =  obj.rhsMaster{i} + ...
        obj.E'*obj.solvers(2).getSolver(obj.physics).rhs;
    end
  end

  methods (Access = public)

    function [dofr,dofc,mat] = computeLocMaster(obj,imult,im,Nmult,Nmaster)
      mat = obj.quadrature.integrate(@(a,b) pagemtimes(a,'ctranspose',b,'none'),...
        Nmult,Nmaster);
      nodeMaster = obj.mesh.local2glob{1}(obj.mesh.msh(1).surfaces(im,:));
      fld = obj.dofm(1).getFieldId(obj.physics);
      dofc = obj.dofm(1).getLocalDoF(nodeMaster,fld);
      dofr = getMultiplierDoF(obj,imult);
    end

    function [dofr,dofc,mat] = computeLocSlave(obj,imult,is,mat)
      if strcmp(obj.multiplierType,'dual')
        % lump local D matrix
        mat = diag(sum(mat,2));
      end
      nodeSlave = obj.mesh.local2glob{2}(obj.mesh.msh(2).surfaces(is,:));
      fld = obj.dofm(2).getFieldId(obj.physics);
      dofc = obj.dofm(2).getLocalDoF(nodeSlave,fld);
      dofr = getMultiplierDoF(obj,imult);
    end


    function [cellStr,pointStr] = buildPrintStruct(obj,fldId,fac)

      fieldName = obj.physics(fldId);
      nCellData = obj.dofm(1).getDoFperEnt(fieldName);
      cellStr = repmat(struct('name', 1, 'data', 1), nCellData, 1);
      if nargin ==2
        outVar = obj.multipliers(fldId).prev;
      elseif nargin == 3
        outVar = fac*obj.multipliers(fldId).curr + ...
          (1-fac)*obj.multipliers(fldId).prev;
      else
        error('Invalid number of input arguments for function buildPrintStruct');
      end

      for i = 1:nCellData
        cellStr(i).name = char(strcat(fieldName,'_',num2str(i)));
        cellStr(i).data = outVar(i:nCellData:end);
      end

      pointStr = [];        % when using P0 multipliers
    end

    function goOnState(obj)
      % update the value of the multipliers
      for i = 1:obj.nFld
        obj.multipliers(i).prev = obj.multipliers(i).curr;
      end
    end

    function goBackState(obj)
      for i = 1:obj.nFld
        obj.multipliers(i).curr = obj.multipliers(i).curr;
      end
    end

    function setDoFcount(obj,ndof)
      obj.dofCount = ndof(obj.idDomain);
    end


    function out = isMatrixComputed(obj)
      out = all(cellfun(@(x) ~isempty(x), [obj.Jmaster(:); obj.Jslave(:)]));
    end

    function Nmult = computeMultiplierBasisF(obj,el,NslaveIn)
      elem = obj.getElem(2,el);
      switch obj.multiplierType
        case 'P0'
          Nmult = ones(size(NslaveIn,1),1);
        case 'standard'
          Nmult = NslaveIn;
        case 'dual'
          Ns = getBasisFinGPoints(elem);
          gpW = getDerBasisFAndDet(elem,el);
          M = Ns'*(Ns.*gpW');
          D = diag(Ns'*gpW');
          A = M\D;
          Nmult = NslaveIn*A;
      end
    end

    function nDofs = getNumbMultipliers(obj)
      nc = obj.dofm(2).getDoFperEnt(obj.physics);
      switch obj.multiplierType
        case 'P0'
          nDofs = nc*obj.mesh.nEl(2);
        case {'dual','standard'}
          nDofs = nc*obj.mesh.msh(2).nNodes;
      end
    end

    function dofs = getMultiplierDoF(obj,id)

      nc = obj.dofm(2).getDoFperEnt(obj.physics);

      if nargin > 1
        if strcmp(obj.multiplierType,'P0')
          dofs = dofId(id,nc);
        else
          is = obj.mesh.activeCells{2}(id);
          nodes = obj.mesh.msh(2).surfaces(is,:);
          dofs = dofId(nodes,nc);
        end
      else
        if strcmp(obj.multiplierType,'P0')
          dofs = dofId(getActiveCells(obj.mesh,2),3);
        else
          dofs = 1:nc*obj.mesh.msh(2).nNodes;
        end
      end
    end

    function [bcDofs,bcVals] = removeSlaveBCdofs(obj,bcPhysics,bcData,domId)
      % this method updates the list of bc dofs and values removing dofs on
      % the slave interface (only for nodal multipliers)
      % this avoid overconstraint and solvability issues

      % bcData: nDofBC x 2 matrix.
      % first column -> dofs
      % second columns -> vals (optional)

      bcDofs = bcData(:,1);
      if size(bcData,2) > 1
        bcVals = bcData(:,2);
      else
        bcVals = [];
      end

      if strcmp(obj.multiplierType,'P0')
        return
      end

      if nargin > 3
        if ~(domId == obj.idDomain(2))
          % not slave side
          return
        end
      end

      % get list of nodal slave dofs in the interface
      nodSlave = obj.mesh.local2glob{2};
      fldId = obj.dofm(2).getFieldId(bcPhysics);
      dofSlave = getLocalDoF(obj.dofm(2),nodSlave,fldId);
      isBCdof = ismember(bcData(:,1),dofSlave);
      bcDofs = bcDofs(~isBCdof);
      if ~isempty(bcVals)
        bcVals = bcVals(~isBCdof);
      end
    end
  end
end

