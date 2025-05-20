classdef MeshGlue < Mortar
  % subclass of Mortar implementing mesh tying between different domains
  % mesh tying enforces continuity between primary variables
  % different physics can be tyied at once

  properties (Access = private)
    mortar
  end

  properties (Access = public)
    physics
    slaveMat
    masterMat
    multipliers
    iniMultipliers
  end

  methods

    function obj = MeshGlue(id,inputStruct,domains)
      obj@Mortar(inputStruct,domains);
      assert(isfield(inputStruct,"Physics"), ...
        'Missing Physics field for interface %i',id);
      obj.physics = split(inputStruct.Physics);
      obj.nFld = numel(obj.physics);
      initializeJacobianAndRhs(obj);
      isFldMaster = isField(obj.dofmMaster,obj.physics);
      isFldSlave = isField(obj.dofmSlave,obj.physics);
      assert(isFldMaster,['MeshGlue physic not available for ' ...
        'master domain %i'],obj.idMaster);
      assert(isFldSlave,['MeshGlue physic not available for ' ...
        'slave domain %i'],obj.idSlave);
      % initializing the multiplier cell array
      computeMortarMatrices(obj);
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
        case 2
          if domId == obj.idMaster
            varargout{1} = (obj.Jmaster{fldId})';
            varargout{2} = obj.Jmaster{fldId};
          elseif domId == obj.idSlave
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
          if domId == obj.idMaster
            rhs = obj.rhsMaster{fldId};
          elseif domId == obj.idSlave
            rhs = obj.rhsSlave{fldId};
          else
            error('Input domain %i is not a valid master/slave',domId)
          end
      end
    end


    function computeMat(obj,idDomain)
      % return matrices for master and slave side in appropriate field
      side = getSide(obj,idDomain);
      for i = 1:obj.nFld
        % map local mortar matrices to global indices
        switch side
          case 'master'
            obj.Jmaster{i} = obj.getMatrix('master',obj.physics{i});
          case 'slave'
            obj.Jslave{i} = obj.getMatrix('slave',obj.physics{i});
        end
      end
    end

    function computeRhs(obj,idDomain,state)
      % compute rhs contributions for a specified input field
      side = getSide(obj,idDomain);
      for i = 1:obj.nFld
        switch side
          case 'master'
            computeRhsMaster(obj,i,state);
          case 'slave'
            computeRhsSlave(obj,i,state);
        end
      end
    end

    function mat = getMatrix(obj,side,field)
      switch side
        case 'master'
          n = obj.dofmMaster.getDoFperEnt(field);
          dofMult = dofId(1:obj.nElSlave,n);
          dofMaster = obj.loc2globMaster(1:size(obj.masterMat,2));
          dofMaster = obj.dofmMaster.getLocalDoF(dofMaster,field);
          [j,i] = meshgrid(dofMaster,dofMult);
          nr = n*obj.nElSlave;
          nc = obj.dofmMaster.getNumDoF(field);
          vals = Discretizer.expandMat(obj.masterMat,n);
          mat = sparse(i(:),j(:),vals(:),nr,nc); % minus sign!
        case 'slave'
          n = obj.dofmSlave.getDoFperEnt(field);
          dofMult = dofId(1:obj.nElSlave,n);
          dofSlave = obj.loc2globSlave(1:size(obj.slaveMat,2));
          dofSlave = obj.dofmSlave.getLocalDoF(dofSlave,field);
          [j,i] = meshgrid(dofSlave,dofMult);
          nr = n*obj.nElSlave;
          nc = obj.dofmSlave.getNumDoF(field);
          vals = Discretizer.expandMat(obj.slaveMat,n);
          mat = sparse(i(:),j(:),vals(:),nr,nc); % minus sign!
      end
    end


    function computeMortarMatrices(obj)
      elemSlave = getElem(obj,'slave');
      [imVec,jmVec,MVec] = allocateMatrix(obj,'master');
      [idVec,jdVec,DVec] = allocateMatrix(obj,'slave');
      % Ns: basis function matrix on slave side
      % Nm: basis function matrix on master side
      % Nmult: basis function matrix for multipliers
      cs = 0; % slave matrix entry counter
      cm = 0; % master matrix entry counter
      
      % and interpolation coordinates)
      for i = 1:obj.nElSlave
        is = obj.activeSlaveCells(i);
        masterElems = find(obj.elemConnectivity(:,is));
        if isempty(masterElems)
          continue
        end
        %Compute Slave quantities
        dJWeighed = elemSlave.getDerBasisFAndDet(is,3);
        posGP = getGPointsLocation(elemSlave,is);
        nSlave = obj.mshIntSlave.surfaces(is,:);
        Nslave = getBasisFinGPoints(elemSlave); % Get slave basis functions
        for im = masterElems'
          nMaster = obj.mshIntMaster.surfaces(im,:);
          [Nm,id] = obj.quadrature.getMasterBasisF(im,posGP); % compute interpolated master basis function
          if any(id)
            % get basis function matrices
            Nm = Nm(id,:);
            Ns = Nslave(id,:);
            Nmult = ones(size(Ns,1),1);
            Nm = Discretizer.reshapeBasisF(Nm,1);
            Ns = Discretizer.reshapeBasisF(Ns,1);
            Nmult = Discretizer.reshapeBasisF(Nmult,1);

            % compute slave mortar matrix
            Dloc = pagemtimes(Nmult,'transpose',Ns,'none');
            [idVec,jdVec,DVec,cs] = Discretizer.computeLocalMatrix( ...
              Dloc,idVec,jdVec,DVec,cs,dJWeighed(id),i,nSlave);

            % compute master mortar matrix
            Mloc = pagemtimes(Nmult,'transpose',Nm,'none');
            [imVec,jmVec,MVec,cm] = Discretizer.computeLocalMatrix( ...
              Mloc,imVec,jmVec,MVec,cm,dJWeighed(id),i,nMaster);

            % sort out gauss points already ised
            dJWeighed = dJWeighed(~id);
            posGP = posGP(~id,:);
            Nslave = Nslave(~id,:);
          end
        end
        if ~all(id)
          % track element not fully projected
          fprintf('GP not sorted for slave elem numb %i \n',is);
          c_ns = c_ns + 1;
        end
      end

      % cut vectors for sparse matrix assembly
      imVec = imVec(1:cm); jmVec = jmVec(1:cm); MVec = MVec(1:cm);
      idVec = idVec(1:cs); jdVec = jdVec(1:cs); DVec = DVec(1:cs);

      % assemble mortar matrices in sparse format
      obj.masterMat = sparse(imVec,jmVec,MVec,...
        obj.nElSlave,obj.mshIntMaster.nNodes);
      obj.slaveMat = sparse(idVec,jdVec,DVec,...
        obj.nElSlave,obj.mshIntSlave.nNodes);
    end

    function applyBCmaster(obj,bound,bc,t,state)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(1).getSolver(physic),bound,bc,t,state);
      i = strcmp(obj.physics,physic);
      obj.rhsMaster{i}(bcEnts) = 0;
      obj.Jmaster{i}(:,bcEnts) = 0;
    end

    function applyBCslave(obj,bound,bc,t,state)
      physic = bound.getPhysics(bc);
      [bcEnts,~] = getBC(obj.solvers(2).getSolver(physic),bound,bc,t,state);
      i = strcmp(obj.physics,physic);
      obj.rhsSlave{i}(bcEnts) = 0;
      obj.Jslave{i}(:,bcEnts) = 0;
    end
  end

  methods (Access = private)

    function initializeJacobianAndRhs(obj)
      [obj.Jmaster,obj.Jslave] = deal(cell(obj.nFld,1));
      [obj.rhsMaster,obj.rhsSlave,obj.multipliers] = deal(cell(obj.nFld,1));
      for i = 1:obj.nFld
        ncomp = obj.dofmMaster.getDoFperEnt(obj.physics(i));
        nDofMaster = getNumDoF(obj.dofmMaster,obj.physics(i));
        nDofSlave = getNumDoF(obj.dofmSlave,obj.physics(i));
        nDofMult = ncomp*obj.nElSlave;
        obj.rhsMaster{i} = zeros(nDofMaster,1);
        obj.rhsSlave{i} = zeros(nDofSlave,1);
        obj.rhsMult{i} = zeros(nDofMult,1);
        obj.multipliers{i} = zeros(nDofMult,1);
        obj.iniMultipliers = obj.multipliers;
      end
    end
    
    function computeRhsMaster(obj,i,state)     
      obj.rhsMaster{i} = ...
        obj.Jmaster{i}'*(obj.multipliers{i}-obj.iniMultipliers{i});
      var = getState(obj.solvers(1).getSolver(obj.physics(i)),state);
      ents = obj.dofmMaster.getActiveEnts(obj.physics(i));
      obj.rhsMult{i} = obj.rhsMult{i} + obj.Jmaster{i}*var(ents);
    end

    function computeRhsSlave(obj,i,state)
      obj.rhsSlave{i} = ...
        obj.Jslave{i}'*(obj.multipliers{i}-obj.iniMultipliers{i});
      var = getState(obj.solvers(2).getSolver(obj.physics(i)),state);
      ents = obj.dofmSlave.getActiveEnts(obj.physics(i));
      obj.rhsMult{i} = obj.rhsMult{i} + obj.Jslave{i}*var(ents);
    end

  end
  
end

