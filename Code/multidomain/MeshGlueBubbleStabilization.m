classdef MeshGlueBubbleStabilization < MeshGlue
  % Mesh glue class implementing pressure jump stabilization (on the slave
  % side)
  % We apply static condensation locally, so this class modifies the inner
  % block of the slave domain

  properties
    JcondSlave  % static condensation block for inner slave domain
    localFaceIndex % slave faces with local index in neighboring cell
    
  end

  methods (Access = public)
    function obj = MeshGlueBubbleStabilization(id,inputStruct,domains)
      obj@MeshGlue(id,inputStruct,domains);
      setFaces(obj);
      computeMortarMatrices(obj);
    end

    function computeMortarMatrices(obj)

      % compute mortar matrices and
      elemSlave = getElem(obj,2);
      [imVec,jmVec,MVec] = allocateMatrix(obj,1);
      [idVec,jdVec,DVec] = allocateMatrix(obj,2);
      % Ns: basis function matrix on slave side
      % Nm: basis function matrix on master side
      % Nmult: basis function matrix for multipliers
      cs = 0; % slave matrix entry counter
      cm = 0; % master matrix entry counter

      for i = 1:obj.mesh.nEl(2)
        is = obj.mesh.activeCells{2}(i);
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

    function computeMat(obj,idDomain)
      computeMat@MeshGlue(obj,idDomain);
    end


    function out = isMatrixComputed(obj)
      out = all(cellfun(@(x) ~isempty(x), ...
        [obj.Jmaster(:); obj.Jslave(:); obj.Jmult(:)]));
    end
  end



  methods(Access = private)

 

    function setFaces(obj)
      % associate the slave element to the local index on the parent cell
      nElSlave = obj.mesh.nEl(2);
      obj.localFaceIndex = zeros(nElSlave,1);
      msh = obj.solvers(2).grid.topology;
      mapf2c = obj.mesh.f2c{2};
      faceNodes = [
        1 2 3 4;
        1 4 5 8;
        1 2 5 6;
        2 3 6 7;
        3 4 7 8;
        5 6 7 8
        ];
      for f = 1:nElSlave
        % get global nodes of the cell
        nodeC = msh.cells(mapf2c(f),:);
        nodeF = obj.mesh.msh(2).surfaces(f,:);
        nodeF = obj.mesh.local2glob{2}(nodeF);
        cellNodes = sort(find(ismember(nodeC,nodeF)));
        % 
        % Face ID lookup using row-wise comparison
        faceID = find(all(bsxfun(@eq, faceNodes, cellNodes), 2), 1);
        obj.localFaceIndex(f) = faceID;
      end
    end




  end
end

