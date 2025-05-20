classdef Mortar < handle
  % This class implement some basic operation for mortar interpolation
  
  properties
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
  end

  methods
    function [obj] = Mortar(inputStruct,domains)
      obj.solvers = [domains(1).Discretizer,domains(2).Discretizer];
      obj.idMaster = inputStruct.Master.idAttribute;
      obj.idSlave = inputStruct.Slave.idAttribute;
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

    function computeStabilizationMatrix(obj,fld)

      % get number of components of input field
      nc = obj.dofmMaster.getDoFperEnt(fld);

      % initialize matrix
      H = zeros(nc*obj.nElSlave);
      
      % get list of internal master edges
      inEdgeMaster = find(all(obj.e2fMaster,2));

      for ie = inEdgeMaster'
        % get edges sharing node n
        edgeID = any(ismember(obj.masterTopol,n),2);
        if sum(edgeID)<2 
          % boundary node
          continue
        else
          eM = find(edgeID); % master cell ID in contact
          % get slave elements in contact
          eS = unique([find(obj.elemConnectivity(eM(1),:)),find(obj.elemConnectivity(eM(2),:))]);
          % node list for these edges
          nSfull = obj.slaveTopol(eS,:);
          [nSglob, ~, ~] = unique(nSfull(:));
          counts = histc(nSfull(:), nSglob);
          % Elements that appear at least two times
          nSglob = nSglob(counts >= 2);
          nSdof = DofMap.getCompDoF(nSglob);
          nS = ismember(obj.nodesSlave,nSglob); % local numbering of slave nodes
          nMfull = obj.masterTopol(eM,:);
          [nMglob, ~, ~] = unique(nMfull(:));
          counts = histc(nMfull(:), nMglob);
          % Elements that appear at least two times
          nMglob = nMglob(counts >= 2);
          nMdof = DofMap.getCompDoF(nMglob);
          nM = ismember(obj.nodesMaster,nMglob); % local numbering of slave nodes
          Km = Kmaster(nMdof,nMdof);
          Ks = Kslave(nSdof,nSdof);
          Dloc = D(eS,nS);
          Mloc = M(eS,nM);
          V = [Dloc'; -Mloc'];
          %Vtilde = zeros(size(V));
          % % swapping columns to get orthogonal matrix
          % Vtilde(:,1:2:end) = -V(:,2:2:end);
          % Vtilde(:,2:2:end) = V(:,1:2:end);
          V = expandMat(V,2);
          Kloc = diag([1./diag(Ks);1./diag(Km)]);
          Hloc = V'*(Kloc*V);
          Hlump = diag(Hloc);
          lM1 = getLength(obj,eM(1),'master');
          lM2 = getLength(obj,eM(2),'master');
          lM = 0.5*(lM1+lM2);
          % assemble stabilization matrix
          for nsLoc = nSglob'
            esLoc = find(any(ismember(obj.slaveTopol(eS,:),nsLoc),2));
            K1 = diag(Hlump([2*esLoc(1)-1 2*esLoc(1)]));
            K2 = diag(Hlump([2*esLoc(2)-1 2*esLoc(2)]));
            K = 0.5*(K1+K2);
            dof1 = [2*esLoc(1)-1 2*esLoc(1)];
            dof2 = [2*esLoc(2)-1 2*esLoc(2)];
            Knew = 0.5*(Hloc(dof1,dof1)+Hloc(dof2,dof2));
            dof = DofMap.getCompDoF(eS(esLoc));
            lS1 = getLength(obj,eS(esLoc(1)),'slave');
            lS2 = getLength(obj,eS(esLoc(2)),'slave');
            lS = 0.5*(lS1+lS2);
            H(dof,dof) = H(dof,dof) + [K -K;-K K];
          end
        end
      end
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
      bot = msh.surfaces(:, [1, 2]);
      top = msh.surfaces(:, [3, 4]);
      lft = msh.surfaces(:, [4, 1]);
      rgt = msh.surfaces(:, [2, 3]);
      % unique matrix containing all existing edges (with repetitions)
      edgeMat = [bot;top;lft;rgt];
      %
      [edgeMat,i] = sort(edgeMat,2);
      i = i(:,1);
      %
      % id is [cellnumber, half-face tag]
      id       = [(1:msh.nSurfaces)', repmat(1, [msh.nSurfaces, 1]);...
        (1:msh.nSurfaces)', repmat(2, [msh.nSurfaces, 1]);...
        (1:msh.nSurfaces)', repmat(3, [msh.nSurfaces, 1]);...
        (1:msh.nSurfaces)', repmat(4, [msh.nSurfaces, 1])]; %#ok<REPMAT>
      % Sort rows to find pairs of cells sharing a face
      [edgeMat, j] = sortrows(edgeMat);

      % encode edge matrix
      [e2n,n] = obj.rle(edgeMat);
      nEdges = numel(n);
      N = repelem(1:nEdges,n);
      switch side
        case 'master'
          obj.e2fMaster = accumarray([N',i(j)],id(j,1),[nEdges,2]);
          obj.e2nMaster = e2n;
          obj.nEdgesMaster = nEdges;
        case 'slave'
          obj.e2fSlave = accumarray([N',i(j)],id(j,1),[nEdges,2]);
          obj.e2nSlave = e2n;
          obj.nEdgesSlave = nEdges;
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

