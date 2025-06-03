classdef MeshGlue < Mortar
  % subclass of Mortar implementing mesh tying between different domains
  % mesh tying enforces continuity between primary variables
  % different physics can be tyied at once

  properties (Access = private)
    mortar
    out
    dofCount
  end

  properties (Access = public)
    physics
    multipliers
    iniMultipliers
    totMult
  end

  methods

    function obj = MeshGlue(id,inputStruct,domains)
      obj@Mortar(inputStruct,domains);
      assert(isfield(inputStruct,"Physics"), ...
        'Missing Physics field for interface %i',id);
      obj.physics = split(inputStruct.Physics);
      obj.nFld = numel(obj.physics);
      initializeJacobianAndRhs(obj);
      isFldMaster = isField(obj.dofm(1),obj.physics);
      isFldSlave = isField(obj.dofm(2),obj.physics);
      assert(isFldMaster,['MeshGlue physic not available for ' ...
        'master domain %i'],obj.idDomain(1));
      assert(isFldSlave,['MeshGlue physic not available for ' ...
        'slave domain %i'],obj.idDomain(2));
      % initializing the multiplier cell array
    end
    %
  end

  methods (Access=public)
    %
    function varargout = getJacobian(obj,fldId,domId)
      % get jacobian blocks associated to specific field and specific
      % domain
      % if nargout = 2 -> get master/slave pair of jacobian blocks
      % if nargout = 1 -> return multiplier jacobian
      switch nargout
        case 1
          varargout{1} = obj.Jmult{fldId};
          % 0.5 is needed because the multiplier matrix is retrieved twice
        case 2
          if domId == obj.idDomain(1)
            varargout{1} = (obj.Jmaster{fldId})';
            varargout{2} = obj.Jmaster{fldId};
          elseif domId == obj.idDomain(2)
            varargout{1} = (obj.Jslave{fldId})';
            varargout{2} = obj.Jslave{fldId};
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end

    function rhs = getRhs(obj,fldId,varargin)
      % return rhs block associated to master/slave/multiplier field
      switch nargin
        case 2
          % multiplier rhs block
          rhs = obj.rhsMult{fldId};
        case 3
          domId = varargin{1};
          if domId == obj.idDomain(1)
            rhs = obj.rhsMaster{fldId};
          elseif domId == obj.idDomain(2)
            rhs = obj.rhsSlave{fldId};
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end


    function computeMat(obj,~)
      % return matrices for master and slave side in appropriate field
      if obj.isMatrixComputed()
        % mesh glue matrices are constant troughout the simulation
        return
      end
      computeMortarMatrices(obj);
      for i = 1:obj.nFld
        % map local mortar matrices to global indices
        obj.Jmaster{i} = obj.getMatrix(1,obj.physics(i));
        obj.Jslave{i} = obj.getMatrix(2,obj.physics(i));
      end
    end

    function computeRhs(obj)
      % compute rhs contributions for a specified input field
      for i = 1:obj.nFld
        computeRhsMaster(obj,i);
        computeRhsSlave(obj,i);
      end
    end

    function applyBCmaster(obj,bound,bc,t,state)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(physic),bound,bc,t);
      i = strcmp(obj.physics,physic);
      obj.rhsMaster{i}(bcEnts) = 0;
      obj.Jmaster{i}(:,bcEnts) = 0;
    end

    function applyBCslave(obj,bound,bc,t,state)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(physic),bound,bc,t);
      i = strcmp(obj.physics,physic);
      obj.rhsSlave{i}(bcEnts) = 0;
      obj.Jslave{i}(:,bcEnts) = 0;
    end

    function updateState(obj,du)
      for i = 1:obj.nFld
        n = numel(obj.multipliers(i).curr);
        obj.multipliers(i).curr = obj.multipliers(i).curr + du(1:n);
        du = du(n+1:end);
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
        ncomp = obj.dofm(1).getDoFperEnt(obj.physics(i));
        nDofMaster = getNumDoF(obj.dofm(1),obj.physics(i));
        nDofSlave = getNumDoF(obj.dofm(2),obj.physics(i));
        nDofMult = ncomp*obj.mesh.nEl(2);
        obj.rhsMaster{i} = zeros(nDofMaster,1);
        obj.rhsSlave{i} = zeros(nDofSlave,1);
        obj.rhsMult{i} = zeros(nDofMult,1);
        obj.multipliers(i).curr = zeros(nDofMult,1);
        obj.multipliers(i).prev = obj.multipliers(i).curr;
        obj.iniMultipliers{i} = obj.multipliers(i).curr;
        obj.totMult = obj.totMult + nDofMult;
      end
    end

    function computeRhsMaster(obj,i)
      obj.rhsMaster{i} = ...
        obj.Jmaster{i}'*(obj.multipliers(i).curr-obj.iniMultipliers{i});
      var = getState(obj.solvers(1).getSolver(obj.physics(i)));
      ents = obj.dofm(1).getActiveEnts(obj.physics(i));
      obj.rhsMult{i} = obj.rhsMult{i} + obj.Jmaster{i}*var(ents);
    end

    function computeRhsSlave(obj,i)
      obj.rhsSlave{i} = ...
        obj.Jslave{i}'*(obj.multipliers(i).curr-obj.iniMultipliers{i});
      var = getState(obj.solvers(2).getSolver(obj.physics(i)));
      ents = obj.dofm(2).getActiveEnts(obj.physics(i));
      obj.rhsMult{i} = obj.rhsMult{i} + obj.Jslave{i}*var(ents); 
      if ~isempty(obj.Jmult{i})
        obj.rhsMult{i} = obj.rhsMult{i} + obj.Jmult{i}*(obj.multipliers(i).curr-obj.iniMultipliers{i});
      end
    end
  end

  methods (Access = public)
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

  end

end

