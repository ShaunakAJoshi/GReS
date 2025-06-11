classdef Mortar < handle
  % This class implement some basic operation for mortar interpolation
  
  properties
    outStruct   
    solvers
    mesh         % instance of interfaceMesh class
    idDomain
    Jmaster
    Jslave
    Jmult
    rhsMaster
    rhsSlave
    rhsMult
    quadrature % integration scheme used for the interface
    % (RBF or ElementBased)
    dofm
    nFld          
    mortarMatrix
    elements
  end

  methods
    function [obj] = Mortar(inputStruct,domains)
      obj.solvers = [domains(1).Discretizer,domains(2).Discretizer];
      obj.idDomain = [inputStruct.Master.idAttribute;
                inputStruct.Slave.idAttribute];
      surfId = {inputStruct.Master.surfaceTagAttribute;
        inputStruct.Slave.surfaceTagAttribute};
      obj.mesh = interfaceMesh(domains,surfId);
      obj.dofm = [domains(1).DoFManager;
                        domains(2).DoFManager];
      switch inputStruct.Quadrature.typeAttribute
        case 'RBF'
          nG = inputStruct.Quadrature.nGPAttribute;
          nInt = inputStruct.Quadrature.nIntAttribute;
          obj.elements = [Elements(obj.mesh.msh(1),1,nG),...
            Elements(obj.mesh.msh(2),1,nG)];
          obj.quadrature = RBFquadrature(obj,nInt);
        case 'ElementBased'
          % Element based will be implemented in the future
        case 'SegmentBased'
          obj.quadrature = SegmentBasedQuadrature(obj,inputStruct.Quadrature.nGPAttribute);
          nG = 2; % dummy nG for elements deifnition
          obj.elements = [Elements(obj.mesh.msh(1),1,nG),...
            Elements(obj.mesh.msh(2),1,nG)];
      end
      setPrintUtils(obj,inputStruct,domains(2).OutState);
    end

    function [r,c,v] = allocateMatrix(obj,sideID)
      nEntries = nnz(obj.mesh.elemConnectivity)...
        *obj.mesh.nN(sideID);
      [r,c,v] = deal(zeros(nEntries,1));
    end


    function applyBC(obj,idDomain,bound,t,state)
      side = getSide(obj,idDomain);
      bcList = bound.db.keys;

      for bc = string(bcList)
        field = bound.getPhysics(bc);
        if ~ismember(field,obj.physics)
          continue
        end

        if ~strcmp(bound.getType(bc),'Dir')
          continue
          % only dirichlet bc has to be enforced to mortar blocks
        else
          switch side
            case 'master'
              applyBCmaster(obj,bound,bc,t,state)
            case 'slave'
              applyBCslave(obj,bound,bc,t,state)
          end
        end
      end
    end

    function elem = getElem(obj,sideID,id)
      % get instance of element class on one the sides of the interface
      % Assumption: same element type in the entire interface
      % get istance of element class based on cell type
      type = obj.mesh.msh(sideID).surfaceVTKType(id);
      elem = obj.elements(sideID).getElement(type);
    end

    function mat = getMatrix(obj,sideID,field)
      n = obj.dofm(sideID).getDoFperEnt(field);
      dofMult = dofId(1:obj.mesh.nEl(2),n);
      dof = obj.mesh.local2glob{sideID}(1:size(obj.mortarMatrix{sideID},2));
      fld = obj.dofm(sideID).getFieldId(field);
      dof = obj.dofm(sideID).getLocalDoF(dof,fld);
      [j,i] = meshgrid(dof,dofMult);
      nr = n*obj.mesh.nEl(2);
      nc = obj.dofm(sideID).getNumDoF(field);
      vals = Discretizer.expandMat(obj.mortarMatrix{sideID},n);
      mat = sparse(i(:),j(:),vals(:),nr,nc); % minus sign!
    end

    function computeMortarMatrices(obj)
      if isempty(obj.mortarMatrix) % compute only once
        elemSlave = getElem(obj,2);
        [imVec,jmVec,MVec] = allocateMatrix(obj,1);
        [idVec,jdVec,DVec] = allocateMatrix(obj,2);
        % Ns: basis function matrix on slave side
        % Nm: basis function matrix on master side
        % Nmult: basis function matrix for multipliers
        cs = 0; % slave matrix entry counter
        cm = 0; % master matrix entry counter

        for i = 1:obj.mesh.nEl(2)
          is = obj.mesh.getActiveCells(2,i);
          masterElems = find(obj.mesh.elemConnectivity(:,is));
          if isempty(masterElems)
            continue
          end
          %Compute Slave quantities
          dJWeighed = elemSlave.getDerBasisFAndDet(is,3);
          posGP = getGPointsLocation(elemSlave,is);
          nSlave = obj.mesh.msh(2).surfaces(is,:);
          Nslave = getBasisFinGPoints(elemSlave); % Get slave basis functions
          for im = masterElems'
            nMaster = obj.mesh.msh(1).surfaces(im,:);
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
            else
              % pair of elements does not share support. update connectivity
              % matrix
              obj.mesh.elemConnectivity(im,is) = 0;
            end
          end
          if ~all(id)
            % track element not fully projected
            fprintf('%i GP not sorted for slave elem numb %i \n',sum(id),is);
          end
        end

        % cut vectors for sparse matrix assembly
        imVec = imVec(1:cm); jmVec = jmVec(1:cm); MVec = MVec(1:cm);
        idVec = idVec(1:cs); jdVec = jdVec(1:cs); DVec = DVec(1:cs);

        % assemble mortar matrices in sparse format
        obj.mortarMatrix{1} = -sparse(imVec,jmVec,MVec,...
          obj.mesh.nEl(2),obj.mesh.msh(1).nNodes);
        obj.mortarMatrix{2} = sparse(idVec,jdVec,DVec,...
          obj.mesh.nEl(2),obj.mesh.msh(2).nNodes);
      end
    end

    function sideStr = getSide(obj,idDomain)
      % get side of the interface 'master' or 'slave' based on the
      % domain input id
      isMaster = obj.idDomain(1) == idDomain;
      isSlave = obj.idDomain(2) == idDomain;
      if isMaster
        sideStr = 'master';
      elseif isSlave
        sideStr = 'slave';
      else
        % consider the case where both sides belong to the same domain,
        % something like 'master_slave'
        error('Input domain not belonging to the interface');
      end
    end

    function finalizeOutput(obj)
      obj.outStruct.VTK.finalize();
    end

    function printState(obj,tOld,tNew)
      cellData2D = [];
      pointData2D = [];
      if nargin == 2
        t = tOld;
        for i = 1:obj.nFld
          [cellData,pointData] = buildPrintStruct(obj,i);
          cellData2D = OutState.mergeOutFields(cellData2D,cellData);
          pointData2D = OutState.mergeOutFields(pointData2D,pointData);
        end
        obj.VTK.writeVTKFile(t, [], [], pointData2D, cellData2D);
      elseif nargin == 3
        tList = obj.outStruct.tList;
        tID = obj.outStruct.tID;
        if tID <= length(tList)
          while tList(tID) <= tNew
            t = tList(tID);
            % Linear interpolation
            fac = (t - tOld)/(tNew - tOld);
            for i = 1:obj.nFld
              [cellData,pointData] = buildPrintStruct(obj,i,fac);
              cellData2D = OutState.mergeOutFields(cellData2D,cellData);
              pointData2D = OutState.mergeOutFields(pointData2D,pointData);
            end
            tID = tID + 1;
            obj.outStruct.VTK.writeVTKFile(t, [], [], pointData2D, cellData2D);
            if tID > length(tList)
              break
            end
          end
          obj.outStruct.tID = tID;
        end
      end
    end

  end

  methods (Access = private)

    function setPrintUtils(obj,str,outState)
      
      if ~isfield(str,'Print')
        return
      else
        out = struct(...
          'name',[],...
          'tID', 1,...
          'tList',[],...
          'VTK',[]);

        out.name = str.Print.nameAttribute;
        out.VTK = VTKOutput(obj.mesh.msh(2),out.name);
        out.tList = outState.timeList;
        obj.outStruct = out;
      end

    end
  end


  methods (Static)
    function [interfaceStruct,modelStruct] = buildInterfaceStruct(fileName,modelStruct)
      fprintf('Mortar initialization... \n')
      % read interface file and construct array of MeshGlue objects
      interfStr = readstruct(fileName);
      interfStr = interfStr.Interface;
      nInterfaces = numel(interfStr);
      interfaceStruct = cell(nInterfaces,1);
      for i = 1:nInterfaces
        idMaster = interfStr(i).Master.idAttribute;
        idSlave = interfStr(i).Slave.idAttribute;
        type = interfStr(i).Type;
        switch type
          case 'MeshTying'
            if strcmp(interfStr.Stabilization,'Jump')
              interfaceStruct{i} = MeshGlueJumpStabilization(i,interfStr(i), ...
                modelStruct([idMaster,idSlave]));
            elseif strcmp(interfStr.Stabilization,'Bubble')
              interfaceStruct{i} = MeshGlueBubbleStabilization(i,interfStr(i), ...
                modelStruct([idMaster,idSlave]));
            end
          case 'Fault'
            % not yet implemented!
            interfaceStruct{i} = Fault();
          otherwise
            error(['Invalid interface law type for interface %i in file' ...
              '%s. \nAvailable types are: \nMeshTying \nFault'],i,fileName);
        end
        addInterface(modelStruct(idMaster).Discretizer,i);
        addInterface(modelStruct(idSlave).Discretizer,i);
      end
      fprintf('Done Mortar initialization. \n')
    end

    function varargout = reshapeBasisFunctions(nc,varargin)
      assert(numel(varargin)==nargout);
      varargout = cell(1,nargout);
      for i = 1:numel(varargin)
        varargout{i} = Mortar.reshapeBasisF(varargin{i},nc);
      end
    end


    function Nout = reshapeBasisF(basis,nc)
      % reshape basis functions to obtain displacement shape function
      % input: nG x nN matrix
      % output: 3 x nN x nG
      [ng,nn,nt] = size(basis);
      Nout = zeros(nc,nc*nn,ng,nt);
      for i = 1:nt
        % index pattern for matrix reshaping
        ind = repmat(linspace(1,nc^2,nc),1,nn)+(nc^2)*repelem(0:nn-1,1,nc);
        N = zeros(nc*nn,ng); % initialize single page
        N(ind(1:nc*nn),:) = repelem(basis(:,:,i)',nc,1);
        Nout(:,:,:,i) = reshape(N,[nc,nn*nc,ng]); % reshaped 3D matrix
      end
    end


  end
end

