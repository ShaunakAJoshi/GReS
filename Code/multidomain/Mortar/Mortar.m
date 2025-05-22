classdef Mortar < handle
  % This class implement some basic operation for mortar interpolation
  
  properties
    name        % identifier for the interface
    printFlag
    outStruct
    solvers
    idMaster
    idSlave
    Jmaster
    Jslave
    Jmult
    rhsMaster
    rhsSlave
    rhsMult
    mshIntMaster
    mshIntSlave
    loc2globMaster
    loc2globSlave
    elemConnectivity
    nGP
    quadrature % integration scheme used for the interface
    % (RBF or ElementBased)
    nNslave
    nNmaster
    nElSlave
    nElMaster
    activeSlaveCells
    activeMasterCells
    dofmMaster
    dofmSlave
    slaveCellType
    masterCellType
    nFld          
    nEdgesMaster
    nEdgesSlave
    e2nMaster     % map edge to connected nodes
    e2nSlave
    e2fMaster     % map edge to connected faces 
    e2fSlave
    f2cMaster     % map face to connected cells
    f2cSlave      
    f2eMaster     % map face to connected edges
    f2eSlave
    mshGlobMaster % entire mesh of master domain 
    mshGlobSlave  % entire mesh of slave domain
    masterMat
    slaveMat
  end

  methods
    function [obj] = Mortar(inputStruct,domains)
      obj.solvers = [domains(1).Discretizer,domains(2).Discretizer];
      obj.idMaster = inputStruct.Master.idAttribute;
      obj.idSlave = inputStruct.Slave.idAttribute;
      obj.name = inputStruct.Name;
      obj.printFlag = inputStruct.Print;
      surfId = {inputStruct.Master.surfaceTagAttribute;
        inputStruct.Slave.surfaceTagAttribute};
      mshMaster = domains(1).Grid.topology;
      mshSlave = domains(2).Grid.topology;
      obj.dofmMaster = domains(1).DoFManager;
      obj.dofmSlave = domains(2).DoFManager;
      obj.mshIntMaster = mshMaster.getSurfaceMesh(surfId{1});
      obj.mshIntSlave = mshSlave.getSurfaceMesh(surfId{2});
      obj.nGP = inputStruct.nGP;
      getCellTypes(obj);
      getConnectivityMatrix(obj);
      mapLocalNod2Glob(obj,'master',mshMaster,surfId{1});
      mapLocalNod2Glob(obj,'slave',mshSlave,surfId{2});
      setupEdgeTopology(obj,'master');
      setupEdgeTopology(obj,'slave');
      obj.buildFace2CellMap(mshMaster,mshSlave);
      obj.mshGlobMaster = mshMaster;
      obj.mshGlobSlave = mshSlave;
      switch inputStruct.Quadrature.typeAttribute
        case 'RBF'
          obj.quadrature = RBF(obj,inputStruct.Quadrature.nIntAttribute);
        case 'ElementBased'
          % Element based will be implemented in the future
      end
    end

    function getConnectivityMatrix(obj)
      cs = ContactSearching(obj.mshIntMaster,obj.mshIntSlave);
      obj.elemConnectivity = cs.elemConnectivity;
      obj.activeSlaveCells = find(any(obj.elemConnectivity,1));
      obj.activeMasterCells = find(any(obj.elemConnectivity,2));
      obj.nElSlave = numel(obj.activeSlaveCells);
      obj.nElMaster = sum(any(obj.elemConnectivity,2)); 
    end

    function [r,c,v] = allocateMatrix(obj,side)
      switch side
        case 'master'
          nEntries = nnz(obj.elemConnectivity)...
            *obj.nNmaster;
        case 'slave'
          nEntries = nnz(obj.elemConnectivity)*...
            obj.nNslave;
      end
      [r,c,v] = deal(zeros(nEntries,1));
    end

    function elem = getElem(obj,side)
      % get instance of element class on one the sides of the interface
      % Assumption: same element type in the entire interface
      switch side
        case 'master'
          % get istance of element class based on cell type
          gM = Gauss(obj.masterCellType,obj.nGP,2);
          elem_tmp = Elements(obj.mshIntMaster,gM);
          if obj.mshIntMaster.surfaceVTKType(1) == 5
            elem = elem_tmp.tri;
          elseif obj.mshIntMaster.surfaceVTKType(1) == 9
            elem = elem_tmp.quad;
          end
        case 'slave'
          gS = Gauss(obj.slaveCellType,obj.nGP,2);
          elem_tmp = Elements(obj.mshIntSlave,gS);
          if obj.mshIntSlave.surfaceVTKType(1) == 5
            elem = elem_tmp.tri;
          elseif obj.mshIntSlave.surfaceVTKType(1) == 9
            elem = elem_tmp.quad;
          end
      end
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
          mat = sparse(i(:),j(:),- vals(:),nr,nc); % minus sign!
        case 'slave'
          n = obj.dofmSlave.getDoFperEnt(field);
          dofMult = dofId(1:obj.nElSlave,n);
          dofSlave = obj.loc2globSlave(1:size(obj.slaveMat,2));
          dofSlave = obj.dofmSlave.getLocalDoF(dofSlave,field);
          [j,i] = meshgrid(dofSlave,dofMult);
          nr = n*obj.nElSlave;
          nc = obj.dofmSlave.getNumDoF(field);
          vals = Discretizer.expandMat(obj.slaveMat,n);
          mat = sparse(i(:),j(:),vals(:),nr,nc); 
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
          else
            % pair of elements does not share support. update connectivity
            % matrix
            obj.elemConnectivity(im,is) = 0;
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

    function stabMat = computeStabilizationMatrix(obj,fld)

      % get number of components of input field
      nc = obj.dofmMaster.getDoFperEnt(fld);

      % initialize matrix estimating number of entries
      % number of internal slave elements
      nes = sum(all(obj.e2fSlave,2));
      nEntries = 2*nc*nes; % each cell should contribute at least two times
      [id1,id2,vals] = deal(zeros(nEntries,1));

      c = 0;

      % get list of internal master edges
      inEdgeMaster = find(all(obj.e2fMaster,2));

      for ieM = inEdgeMaster'
        % get master faces sharing internal edge ie
        fM = obj.e2fMaster(ieM,:);
        assert(numel(fM)==2,['Unexpected number of connected faces for' ...
          'master edge %i. Expected 2.'], ieM);

        % get slave faces sharing support with master faces
        fS = unique([find(obj.elemConnectivity(fM(1),:)),...
          find(obj.elemConnectivity(fM(2),:))]);

        % get internal edges of slave faces
        eS = unique(obj.f2eSlave(fS,:));
        id = all(ismember(obj.e2fSlave(eS,:),fS),2);
        ieS = eS(id);

        % get active macroelement nodes
        nM = obj.e2nMaster(ieM,:);
        nS = unique(obj.e2nSlave(eS,:));

        % compute local schur complement approximation
        S = computeSchurLocal(obj,nM,nS,fS,fld);

        % assemble stabilization matrix component
        for iesLoc = ieS'
           f = obj.e2fSlave(iesLoc,:);
           fL = find(fS==f(1)); fR = find(fS==f(2));
           Sdiag = diag(S);
           K = 0.5*(Sdiag(dofId(fL,nc))+Sdiag(dofId(fR,nc)));
           id1(c+1:c+nc) = dofId(f(1),nc);
           id2(c+1:c+nc) = dofId(f(2),nc);
           vals(c+1:c+nc) = K;
           c = c+nc;
        end
      end
      
      id1 = id1(1:c); id2 = id2(1:c); vals = vals(1:c);
      % assemble sparse matrix
      nmult = nc*obj.nElSlave;
      stabMat = sparse(id1,id1,vals,nmult,nmult)+...
        sparse(id1,id2,-vals,nmult,nmult)+...
        sparse(id2,id2,vals,nmult,nmult);
      stabMat = stabMat + stabMat' - diag(diag(stabMat));
    end

    function S = computeSchurLocal(obj,nm,ns,fs,fld)
      % compute approximate schur complement for local nonconforming
      % patch of element
      % input: nm/ns local master/slave node indices
      % fs: local slave faces indices
      
      nc = obj.dofmMaster.getDoFperEnt(fld);

      % get local mortar matrices
      Dloc = obj.slaveMat(fs,ns);
      Mloc = obj.masterMat(fs,nm);
      V = [Dloc, -Mloc];              % minus sign!
      V = Discretizer.expandMat(V,nc);

      % get slave and master dof to access jacobian
      dofS = obj.dofmSlave.getLocalDoF(obj.loc2globSlave(ns),fld);
      dofM = obj.dofmMaster.getLocalDoF(obj.loc2globMaster(nm),fld);

      % get local jacobian
      Km = getSolver(obj.solvers(1),fld).J(dofM,dofM);
      Ks = getSolver(obj.solvers(2),fld).J(dofS,dofS);
      Kloc = diag([1./diag(Ks);1./diag(Km)]);

      S = V*(Kloc*V');  % compute Schur complement
    end

    function sideStr = getSide(obj,idDomain)
      % get side of the interface 'master' or 'slave' based on the
      % domain input id
      isMaster = obj.idMaster == idDomain;
      isSlave = obj.idSlave == idDomain;
      if isMaster
        sideStr = 'master';
      elseif isSlave
        sideStr = 'slave';
      else
        error('Input domain not belonging to the interface');
      end
    end
  end

  methods (Access = private)

    function setupEdgeTopology(obj,side)
      % reorder surface topology
      % inspired by face topology for FV in MRST

      switch side
        case 'master'
          msh = obj.mshIntMaster;
        case 'slave'
          msh = obj.mshIntSlave;
      end

      % Extract edges from each face (assuming quads)
      bot = msh.surfaces(:, [1, 2]);  % bottom edge
      top = msh.surfaces(:, [3, 4]);  % top edge
      lft = msh.surfaces(:, [4, 1]);  % left edge
      rgt = msh.surfaces(:, [2, 3]);  % right edge

      % Stack all edges (with duplicates)
      edgeMat = [bot; top; lft; rgt];

      % Sort nodes within each edge to ensure consistency (edge = {min,max})
      [edgeMat, i] = sort(edgeMat, 2);  % sort each row
      i = i(:,1);  % keep first index (for left/right choice later)

      % Build [faceID, localEdgeIndex] for each edge occurrence
      nF = msh.nSurfaces;
      id = [ (1:nF)', repmat(1, nF, 1);    % bottom
        (1:nF)', repmat(2, nF, 1);    % top
        (1:nF)', repmat(3, nF, 1);    % left
        (1:nF)', repmat(4, nF, 1) ];  % right

      % Sort rows of edgeMat to group identical edges
      [edgeMat, j] = sortrows(edgeMat);

      % Run-length encode to find unique edges
      [e2n, n] = obj.rle(edgeMat);
      nEdges = numel(n);

      % For each occurrence, get the corresponding edge index
      N = repelem(1:nEdges, n);  % same size as edgeMat

      % Build face-to-edge mapping: each face maps to 4 edge indices
      edgeIDs = zeros(size(edgeMat, 1), 1);
      edgeIDs(j) = N;  % reorder to match original edge order

      % Now reshape edgeIDs into f2e: one row per face, 4 edges per face
      f2e = reshape(edgeIDs, [nF, 4]);  % columns: [bot, top, left, right]

      % Assign fields depending on side
      switch side
        case 'master'
          obj.e2fMaster = accumarray([N', i(j)], id(j,1), [nEdges, 2]);
          obj.e2nMaster = e2n;
          obj.nEdgesMaster = nEdges;
          obj.f2eMaster = f2e;
        case 'slave'
          obj.e2fSlave = accumarray([N', i(j)], id(j,1), [nEdges, 2]);
          obj.e2nSlave = e2n;
          obj.nEdgesSlave = nEdges;
          obj.f2eSlave = f2e;
      end
    end

    %
    function getCellTypes(obj)
      switch obj.mshIntSlave.surfaceVTKType(1)
        case 5 % Triangle mesh Master
          obj.slaveCellType = 10;
          obj.nNslave = 3;
        case 9 % Quad mesh Master
          obj.slaveCellType = 12;
          obj.nNslave = 4;
      end

      switch obj.mshIntMaster.surfaceVTKType(1)
        case 5 % Triangle mesh Master
          obj.nNmaster = 3;
          obj.masterCellType = 10;
        case 9 % Quad mesh Master
          obj.nNmaster = 4;
          obj.masterCellType = 12;
      end
    end


    function mapLocalNod2Glob(obj,side,msh,surf)
      % return array where local node numbering map to node id in the full
      % 3D grid
      surfGlob2loc = find(ismember(msh.surfaceTag,surf));
      globNodes = (msh.surfaces(surfGlob2loc,:))';
      globNodes = globNodes(:);
      switch side
        case 'master'
          obj.loc2globMaster = zeros(obj.mshIntMaster.nNodes,1);
          locNodes = (obj.mshIntMaster.surfaces)';
          locNodes = locNodes(:);
          obj.loc2globMaster(locNodes) = globNodes;
        case 'slave'
          obj.loc2globSlave = zeros(obj.mshIntSlave.nNodes,1);
          locNodes = (obj.mshIntSlave.surfaces)';
          locNodes = locNodes(:);
          obj.loc2globSlave(locNodes) = globNodes;
      end
    end

    function buildFace2CellMap(obj,mshM,mshS)
      % for now, this method only work for hexahedral meshes
      slaveTop = obj.loc2globSlave(obj.mshIntSlave.surfaces);
      masterTop = obj.loc2globMaster(obj.mshIntMaster.surfaces);
      % get cell ID containing a face
      obj.f2cMaster = zeros(obj.nElMaster,1);
      obj.f2cSlave = zeros(obj.nElSlave,1);
      listCell = (1:mshS.nCells)';
      for i = 1:obj.nElSlave
        idMat = ismember(mshS.cells(listCell,:),slaveTop(i,:));
        idMat = sum(idMat,2);
        id = find(idMat==4);
        assert(isscalar(id),'Invalid connectivity between face and cells');
        obj.f2cSlave(i) = listCell(id);
        listCell(id) = [];
      end
      listCell = (1:mshM.nCells)';
      for i = 1:obj.nElMaster
        idMat = ismember(mshM.cells(listCell,:),masterTop(i,:));
        idMat = sum(idMat,2);
        id = find(idMat==4);
        assert(isscalar(id),'Invalid connectivity between face and cells');
        obj.f2cMaster(i) = listCell(id);
        listCell(id) = [];
      end
    end
  end


  methods (Static)
    function [interfaceStruct,modelStruct] = buildInterfaceStruct(fileName,modelStruct)

      % read interface file and construct array of MeshGlue objects
      interfStr = readstruct(fileName);
      interfStr = interfStr.Interface;
      nInterfaces = numel(interfStr);
      interfaceStruct = cell(nInterfaces,1);
      for i = 1:nInterfaces
        idMaster = interfStr(i).Master.idAttribute;
        idSlave = interfStr(i).Slave.idAttribute;
        type = interfStr(i).Type;
        if strcmp(type,'MeshTying')
          interfaceStruct{i} = MeshGlue(i,interfStr(i), ...
            modelStruct([idMaster,idSlave]));
        elseif strcmp(type,'Fault')
          % not yet implemented!
          interfaceStruct{i} = Fault();
        else
          error(['Invalid interface law type for interface %i in file' ...
            '%s. \nAvailable types are: \nMeshTying \nFault'],i,fileName);
        end
        addInterface(modelStruct(idMaster).Discretizer,i);
        addInterface(modelStruct(idSlave).Discretizer,i);
      end
    end

    function [outMat,count] = rle(inMat)
      % the input matrix has already been sorted by rows
      % Find unique rows and their first occurrences
      [outMat, firstIdx, ~] = unique(inMat, 'rows', 'stable');

      % Compute counts of each unique row
      count = diff([firstIdx; size(inMat,1) + 1]);

      % Restore original order of indices
      %indices = sortIdx(firstIdx);
    end
  end
end

