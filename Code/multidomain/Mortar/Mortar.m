classdef Mortar < handle
  % This class implement some basic operation for mortar interpolation
  
  properties
    mshMaster
    mshSlave
    slaveIdLoc
    slaveIdGlob
    polytopSize = 18 % standard 3D polytop
    elemConnectivity
    nGP
    quadrature % integration scheme used for the interface
    % (RBF or ElementBased)
    nNslave
    nNmaster
    nElSlave
    nElMaster
    slaveCellType
    masterCellType
  end

  methods
    function obj = Mortar(mG,nG,quadScheme)
      obj.mshMaster = mG.meshMaster.getSurfaceMesh(mG.surfId{1});
      obj.mshSlave = mG.meshSlave.getSurfaceMesh(mG.surfId{2});
      obj.nGP = nG;
      getCellTypes(obj);
      getConnectivityMatrix(obj,mG);
      switch quadScheme.typeAttribute
        case 'RBF'
          obj.quadrature = RBF(obj,quadScheme.nIntAttribute);
        case 'ElementBased'
          % Element based will be implemented in the future
      end
    end

    function getConnectivityMatrix(obj,mG)
      cs = ContactSearching(obj.mshMaster,obj.mshSlave,obj.polytopSize);
      obj.nElMaster = sum(any(cs.elemConnectivity, 2));
      obj.nElSlave = sum(any(cs.elemConnectivity, 1));
      obj.elemConnectivity = cs.elemConnectivity;
      [idM,idS,~] = find(cs.elemConnectivity);
      obj.slaveIdLoc = idS([true; diff(idS) ~= 0]); % faster than unique(idS)
      obj.masterIdLoc = unique(idM);
      % map indices of local surface id to id in global mesh
      idSurfMaster = find(ismember(mG.meshMaster.surfaceTag,mG.surfId{1}));
      idSurfSlave = find(ismember(mG.meshSlave.surfaceTag,mG.surfId{2}));
      obj.slaveIdGlob = idSurfSlave(obj.slaveIdLoc);
      obj.masterIdGlob = idSurfMaster(obj.masterIdLoc);
    end

    function [r,c,v] = allocateMatrix(obj,side,nc)
      switch side
        case 'master'
          nEntries = (nc^2)*nnz(obj.elemConnectivity)...
            *obj.nNmaster^2;
        case 'slave'
          nEntries = (nc^2)*nnz(obj.elemConnectivity)*...
            obj.nNslave^2;
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
          elem_tmp = Elements(obj.mshMaster,gM);
          if obj.mshMaster.surfaceVTKType(1) == 5
            elem = elem_tmp.tri;
          elseif obj.mshMaster.surfaceVTKType(1) == 9
            elem = elem_tmp.quad;
          end
        case 'slave'
          gS = Gauss(obj.slaveCellType,obj.nGP,2);
          elem_tmp = Elements(obj.mshSlave,gS);
          if obj.mshSlave.surfaceVTKType(1) == 5
            elem = elem_tmp.tri;
          elseif obj.mshSlave.surfaceVTKType(1) == 9
            elem = elem_tmp.quad;
          end
      end
    end

  end

  methods (Access = private)
    %
    function getCellTypes(obj)
      switch obj.mshSlave.surfaceVTKType(1)
        case 5 % Triangle mesh Master
          obj.slaveCellType = 10;
          obj.nNslave = 3;
        case 9 % Quad mesh Master
          obj.slaveCellType = 12;
          obj.nNslave = 4;
      end

      switch obj.mshMaster.surfaceVTKType(1)
        case 5 % Triangle mesh Master
          obj.nNmaster = 3;
          obj.masterCellType = 10;
        case 9 % Quad mesh Master
          obj.nNmaster = 4;
          obj.masterCellType = 12;
      end
    end
  end
end

