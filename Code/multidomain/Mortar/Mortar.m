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
    entityMapMaster
    entityMapSlave
    elemConnectivity
    nGP
    quadrature % integration scheme used for the interface
    % (RBF or ElementBased)
    nNslave
    nNmaster
    nElSlave
    nElMaster
    dofmMaster
    dofmSlave
    slaveCellType
    masterCellType
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
      obj.entityMapMaster = find(ismember(mshMaster.surfaceTag,surfId{1}));
      obj.entityMapSlave = find(ismember(mshSlave.surfaceTag,surfId{2}));

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
      [obj.nElMaster, obj.nElSlave] = size(obj.elemConnectivity);
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
  end
end

