classdef MeshGlue < handle
  %
  properties (Access = private)
    mortar
  end

  properties (Access = public)
    domains % 1->master master 2->slave
    idDomains
    meshMaster
    meshSlave
    surfId
    physic
    nComp
  end

  methods
    %
    function setMeshGlue(obj,inputStruct,domains)
      % domains(1): master domains(2): slave
      obj.domains = domains;
      obj.idDomains = [inputStruct.Master.idAttribute;
        inputStruct.Slave.idAttribute];
      obj.surfId = {inputStruct.Master.surfaceTagAttribute;
        inputStruct.Slave.surfaceTagAttribute};
      obj.meshMaster = domains(1).Grid.topology;
      obj.meshSlave = domains(2).Grid.topology;
      obj.physic = inputStruct.Physics;
      obj.nComp = domains(1).DoFManager.getDoFperEnt(obj.physic);
      setMortar(obj,inputStruct);
      %[D,M] = computeMortarMatrices(obj);
    end
  end

  methods (Access=public)
    %
    function [D,M] = computeMortarMatrices(obj)
      elc = 0;
      nDofSlave = obj.nComp*numel(obj.mortar.slaveIdGlob);
      elemSlave = getElem(obj.mortar,'slave');
      [imVec,jmVec,MVec] = allocateMatrix(obj.mortar,'master',obj.nComp);
      [idVec,jdVec,DVec] = allocateMatrix(obj.mortar,'slave',obj.nComp);
      % Ns: basis function matrix on slave side
      % Nm: basis function matrix on master side
      % Nmult: basis function matrix for multipliers
      cs = 0; % slave matrix entry counter
      cm = 0; % master matrix entry counter

      % and interpolation coordinates)
      for is = 1:obj.mortar.nElSlave
        %Compute Slave quantities
        isLoc = obj.mortar.slaveIdLoc(is);
        isGlob = obj.mortar.slaveIdGlob(is);
        dJWeighed = elemSlave.getDerBasisFAndDet(isLoc,3);
        posGP = getGPointsLocation(elemSlave,isLoc);
        nSlave = obj.meshSlave.surfaces(isGlob,:);
        Nslave = getBasisFinGPoints(elemSlave); % Get slave basis functions
        masterElems = find(obj.mortar.elemConnectivity(:,isLoc));
        for im = masterElems'
          nMaster = obj.meshMaster.surfaces(im,:);
          [Nm,id] = obj.mortar.quadrature.getMasterBasisF(im,posGP); % compute interpolated master basis function
          if any(id)
            % get basis function matrix
            Nm = Nm(id,:);
            Ns = Nslave(id,:);
            Nmult = ones(size(Ns,1),1);
            Nm = Discretizer.reshapeBasisF(Nm,obj.nComp);
            Ns = Discretizer.reshapeBasisF(Ns,obj.nComp);
            Nmult = Discretizer.reshapeBasisF(Nmult,obj.nComp);

            % get degrees of freedom of local matrix
            dofMaster = getDof(obj,nMaster,'master');
            dofSlave = getDof(obj,nSlave,'slave');
            dofMult = dofId(is,obj.nComp);

            % compute slave mortar matrix
            Dloc = pagemtimes(Nmult,'transpose',Ns,'none');
            [idVec,jdVec,DVec,cs] = Discretizer.computeLocalMatrix( ...
              Dloc,idVec,jdVec,DVec,cs,dJWeighed(id),dofMult,dofSlave);

            % compute master mortar matrix
            Mloc = pagemtimes(Nmult,'transpose',Nm,'none');
            [imVec,jmVec,MVec,cm] = Discretizer.computeLocalMatrix( ...
              Mloc,imVec,jmVec,MVec,cm,dJWeighed(id),dofMult,dofMaster);

            % sort out gauss points already ised
            dJWeighed = dJWeighed(~id);
            posGP = posGP(~id,:);
            Nslave = Nslave(~id,:);
          else
            elc = elc+1;
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
      M = sparse(imVec,jmVec,MVec,...
        nDofSlave,obj.domains(1).DoFManager.totDoF);
      D = sparse(idVec,jdVec,DVec,...
        nDofSlave,obj.domains(2).DoFManager.totDoF);
    end

    function setMortar(obj,input)
      obj.mortar = Mortar(obj,input.nGP,input.Quadrature);
    end


    function dofs = getDof(obj,list,side)
      % get dof corresponding to entities of
      switch side
        case 'master'
          dofs = getDoF(obj.domains(1).DoFManager,obj.physic,list);
        case 'slave'
          dofs = getDoF(obj.domains(2).DoFManager,obj.physic,list);
      end
    end
  end


  methods (Access = protected)

  end


  methods (Static)
    function interfaceStruct = buildInterfaceStruct(fileName,modelStruct)
      % read interface file and construct array of MeshGlue objects
      interfStr = readstruct(fileName);
      interfStr = interfStr.Interface;
      nInterfaces = numel(interfStr);
      interfaceStruct(nInterfaces, 1) = MeshGlue();  % Preallocate array
      for i = 1:nInterfaces
        idMaster = interfStr(i).Master.idAttribute;
        idSlave = interfStr(i).Slave.idAttribute;
        type = interfStr(i).Type;
        if strcmp(type,'MeshTying')
          interfaceStruct(i).setMeshGlue(interfStr(i), ...
            modelStruct([idMaster,idSlave]));
        elseif strcmp(type,'Fault')
          % not yet implemented!
          interfaceStruct(i) = Fault();
          interfaceStruct(i).setFault(interfStr(i), ...
            modelStruct([idMaster,idSlave]));
        else
          error(['Invalid interface law type for interface %i in file' ...
            '%s. \nAvailable types are: \nMeshTying \nFault'],i,fileName);
        end
      end
    end
  end
end

