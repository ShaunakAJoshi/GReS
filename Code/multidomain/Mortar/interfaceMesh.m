classdef interfaceMesh < handle
  % Manager for topological information and operation on interfaces between
  % meshes

  
  properties
    % prop{1} -> master     prop{2} -> slave 
    msh
    local2glob
    elemConnectivity
    nN
    nEl
    activeCells
    cellType
    nEdges
    e2n
    e2f
    f2c
    f2e
  end
  
  methods
    function obj = interfaceMesh(domains,surf)

      mshMaster = domains(1).Grid.topology;
      mshSlave = domains(2).Grid.topology;
      obj.msh = [mshMaster.getSurfaceMesh(surf{1});
                 mshSlave.getSurfaceMesh(surf{2})];
      getCellTypes(obj);
      getConnectivityMatrix(obj);
      mapLocalNod2Glob(obj,1,mshMaster,surf{1});
      mapLocalNod2Glob(obj,2,mshSlave,surf{2});
      setupEdgeTopology(obj,1);
      setupEdgeTopology(obj,2);
      obj.buildFace2CellMap([mshMaster mshSlave]);
    end

    function getConnectivityMatrix(obj)
      cs = ContactSearching(obj.msh(1),obj.msh(2));
      obj.elemConnectivity = cs.elemConnectivity;
      obj.activeCells{2} = find(any(obj.elemConnectivity,1));
      obj.activeCells{1} = find(any(obj.elemConnectivity,2));
      obj.nEl(2) = numel(obj.activeCells{2});
      obj.nEl(1) = sum(any(obj.elemConnectivity,2));
    end
  end


    methods (Access = private)

    function setupEdgeTopology(obj,sideID)
      % reorder surface topology
      % inspired by face topology for FV in MRST

      mesh = obj.msh(sideID);

      % Extract edges from each face (assuming quads)
      bot = mesh.surfaces(:, [1, 2]);  % bottom edge
      top = mesh.surfaces(:, [3, 4]);  % top edge
      lft = mesh.surfaces(:, [4, 1]);  % left edge
      rgt = mesh.surfaces(:, [2, 3]);  % right edge

      % Stack all edges (with duplicates)
      edgeMat = [bot; top; lft; rgt];

      % Sort nodes within each edge to ensure consistency (edge = {min,max})
      [edgeMat, i] = sort(edgeMat, 2);  % sort each row
      i = i(:,1);  % keep first index (for left/right choice later)

      % Build [faceID, localEdgeIndex] for each edge occurrence
      nF = mesh.nSurfaces;
      id = [ (1:nF)', repmat(1, nF, 1);    % bottom
        (1:nF)', repmat(2, nF, 1);    % top
        (1:nF)', repmat(3, nF, 1);    % left
        (1:nF)', repmat(4, nF, 1) ];  % right

      % Sort rows of edgeMat to group identical edges
      [edgeMat, j] = sortrows(edgeMat);

      % Run-length encode to find unique edges
      [obj.e2n{sideID}, n] = rle(edgeMat);
      obj.nEdges(sideID) = numel(n);

      % For each occurrence, get the corresponding edge index
      N = repelem(1:obj.nEdges(sideID), n);  % same size as edgeMat

      % Build face-to-edge mapping: each face maps to 4 edge indices
      edgeIDs = zeros(size(edgeMat, 1), 1);
      edgeIDs(j) = N;  % reorder to match original edge order

      % Now reshape edgeIDs into f2e: one row per face, 4 edges per face
      obj.f2e{sideID} = reshape(edgeIDs, [nF, 4]);  % columns: [bot, top, left, right]
      
      obj.e2f{sideID} = accumarray([N', i(j)], id(j,1), [obj.nEdges(sideID), 2]);
    end

    %
    function getCellTypes(obj)
      for i = [1 2]
        switch obj.msh(i).surfaceVTKType(1)
          case 5 % Triangle mesh Master
            obj.cellType(i) = 10;
            obj.nN(i) = 3;
          case 9 % Quad mesh Master
            obj.cellType(i) = 12;
            obj.nN(i) = 4;
        end
      end
    end


    function mapLocalNod2Glob(obj,side,msh,surf)
      % return array where local node numbering map to node id in the full
      % 3D grid
      surfGlob2loc = find(ismember(msh.surfaceTag,surf));
      globNodes = (msh.surfaces(surfGlob2loc,:))';
      globNodes = globNodes(:);
      obj.local2glob{side} = zeros(obj.msh(side).nNodes,1);
      locNodes = (obj.msh(side).surfaces)';
      locNodes = locNodes(:);
      obj.local2glob{side}(locNodes) = globNodes;
    end

    function buildFace2CellMap(obj,meshBg)
      % for now, this method only work for hexahedral meshes
      for i = [1 2]
      top = obj.local2glob{i}(obj.msh(i).surfaces);
      % get cell ID containing a face
      obj.f2c{i} = zeros(obj.nEl(i),1);
      listCell = (1:meshBg(i).nCells)';
      for e = 1:obj.nEl(i)
        idMat = ismember(meshBg(i).cells(listCell,:),top(e,:));
        idMat = sum(idMat,2);
        id = find(idMat==4);
        assert(isscalar(id),'Invalid connectivity between face and cells');
        obj.f2c{i}(e) = listCell(id);
        listCell(id) = [];
      end
      end
    end
  end
end

