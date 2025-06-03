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
      assert(obj.nFld==1 && strcmp(obj.physics,'Poromechanics'),['' ...
        'Bubble stabilization is currently available for Poromechanics only'])
      setFaces(obj);
    end

    function computeMortarMatrices(obj,dt)
      % assumption: each element has only one bubble face


      % loop over slave faces and:
      % 1) compute Aub, Abu and Abb local matrix from the neighbor cell
      % 2) compute local M, D and Db
      % 3) assemble static condensation blocks and Jmult 

      % differently from the base method of the mortar class, here global
      % dof indexing is used for the assembled matrices

      % get element instance for mortar integration
      elemSlave = getElem(obj,2);

      % allocate mortar matrices
      Nm = nnz(obj.mesh.elemConnectivity)*9*obj.mesh.nN(1);
      [iMvec,jMvec,Mvec] = deal(zeros(Nm,1));
      %[iDvec,jDvec,Dvec] = allocateMatrix(obj,2);

      Nd = obj.mesh.nEl(2)*9*obj.mesh.nN(2);
      [iDvec,jDvec,Dvec] = deal(zeros(Nd,1));

      fld = [obj.dofm(1).getFieldId(obj.physics), ...
             obj.dofm(2).getFieldId(obj.physics)];
      
      % allocate condensation block for the inner slave domain
      poro = getSolver(obj.solvers(2),obj.physics);
      subCells = obj.mesh.f2c{2};
      mshSlave = obj.solvers(2).grid.topology;
      nSubCellsByType = histc(mshSlave.cellVTKType(subCells),[10, 12, 13, 14]);
      nuu = poro.nEntryKLoc*nSubCellsByType;
      [iKvec,jKvec,Kvec] = deal(zeros(nuu,1));

      % allocate condensation block for coupling slave dof and multipliers
      nul = 9*8*obj.mesh.nEl(2); 
      [iBvec,jBvec,Bvec] = deal(zeros(nul,1));

      % allocate condensation block for multipliers
      [iLvec,jLvec,Lvec] = deal(zeros(9*obj.mesh.nEl(2),1));

      % Ns: basis function matrix on slave side
      % Nm: basis function matrix on master side
      % Nmult: basis function matrix for multipliers
      cd = 0; % slave matrix entry counter
      cm = 0; % master matrix entry counter
      cl = 0; % multiplier matrix counter
      ck = 0; % slave inner block entry counter
      cb = 0; % condensation disp-mult entry counter

      % set number of dofs for each block
      nDofMaster = obj.dofm(1).getNumDoF(obj.physics);
      nDofSlave = obj.dofm(2).getNumDoF(obj.physics);
      nDofMult = obj.mesh.nEl(2)*obj.dofm(2).getDoFperEnt(obj.physics);


      for i = 1:obj.mesh.nEl(2)
        is = obj.mesh.activeCells{2}(i);
        masterElems = find(obj.mesh.elemConnectivity(:,is));
        if isempty(masterElems)
          continue
        end
        %Compute Mortar Slave quantities
        w = elemSlave.getDerBasisFAndDet(is,3);
        posGP = getGPointsLocation(elemSlave,is);
        nodeSlave = obj.mesh.local2glob{2}(obj.mesh.msh(2).surfaces(is,:));
        Nslave = getBasisFinGPoints(elemSlave); % Get slave basis functions
        Nbubble = getBubbleBasisFinGPoints(elemSlave); 
        Dbloc = zeros(3,3);
        Dloc = zeros(3,3*size(Nslave,2));
        for im = masterElems'
          nodeMaster = obj.mesh.local2glob{1}(obj.mesh.msh(1).surfaces(im,:));
          [Nm,id] = obj.quadrature.getMasterBasisF(im,posGP); % compute interpolated master basis function
          if any(id)
            % get basis function matrices
            Nm = Nm(id,:);
            Ns = Nslave(id,:);
            Nb = Nbubble(id,:);
            Nmult = ones(size(Ns,1),1);

            % assuming poromechanics only
            Nm = Discretizer.reshapeBasisF(Nm,3);
            Ns = Discretizer.reshapeBasisF(Ns,3);
            Nb = Discretizer.reshapeBasisF(Nb,3);
            Nmult = Discretizer.reshapeBasisF(Nmult,3);

            % compute slave mortar matrix
            Dtmp = pagemtimes(Nmult,'transpose',Ns,'none');
            Dtmp = Dtmp.*reshape(w(id),1,1,[]);
            Dtmp = sum(Dtmp,3);
            Dloc = Dloc+Dtmp;

            % compute master mortar matrix
            Mtmp = pagemtimes(Nmult,'transpose',Nm,'none');
            Mtmp = Mtmp.*reshape(w(id),1,1,[]);
            Mloc = sum(Mtmp,3);
            nm = numel(Mloc);
            dofMaster = obj.dofm(1).getLocalDoF(nodeMaster,fld(1));
            [jMloc, iMloc] = meshgrid(dofMaster, dofId(i,3));
            iMvec(cm+1:cm+nm) = iMloc(:);
            jMvec(cm+1:cm+nm) = jMloc(:);
            Mvec(cm+1:cm+nm) = Mloc(:);
            cm = cm+nm;
            % [iMvec,jMvec,Mvec,cm] = Discretizer.computeLocalMatrix( ...
            % Mloc,iMvec,jMvec,Mvec,cm,w(id),i,nMaster);

            % compute slave Bubble matrix
            Dbtmp = pagemtimes(Nmult,'transpose',Nb,'none');
            Dbtmp = Dbtmp.*reshape(w(id),1,1,[]);
            Dbtmp = sum(Dbtmp,3);
            Dbloc = Dbloc+Dbtmp;

            % sort out gauss points already used
            w = w(~id);
            posGP = posGP(~id,:);
            Nslave = Nslave(~id,:);
          else
            % pair of elements does not share support. update connectivity
            % matrix
            obj.mesh.elemConnectivity(im,is) = 0;
          end
          if all(id)
            break
          end
        end

        % slave inner condensation contribution
        cellId = obj.mesh.f2c{2}(is);
        [dofRow,dofCol,Kub,Kbb] = computeLocalStiffBubble(poro,cellId,dt);
        [jKloc,iKloc] = meshgrid(dofCol,dofRow);
        % extract only active bubble dofs
        faceId = dofId(obj.localFaceIndex(i),3);
        Kub = Kub(:,faceId);
        Kbb = Kbb(faceId,faceId);
        invKbb = inv(Kbb);
        KuuLoc = -Kub*(invKbb)*Kub';
        nk = numel(KuuLoc);
        iKvec(ck+1:ck+nk) = iKloc(:);   jKvec(ck+1:ck+nk) = jKloc(:);
        Kvec(ck+1:ck+nk) = KuuLoc(:);
        
        % condensation term coupling displacement with multipliers
        Bloc = -Dbloc*(invKbb)*Kub';
        [jBloc,iBloc] = meshgrid(dofRow,dofId(i,3));
        nb = numel(Bloc);
        iBvec(cb+1:cb+nb) = iBloc(:);   jBvec(cb+1:cb+nb) = jBloc(:);
        Bvec(cb+1:cb+nb) = Bloc(:);

        % Mortar Dslave matrix
        dofSlave = obj.dofm(2).getLocalDoF(nodeSlave,fld(2));
        [jDloc,iDloc] = meshgrid(dofSlave,dofId(i,3));
        nd = numel(Dloc);
        iDvec(cd+1:cd+nd) = iDloc(:);   jDvec(cd+1:cd+nd) = jDloc(:);
        Dvec(cd+1:cd+nd) = Dloc(:);

        % multipliers condensation
        Lloc = - Dbloc*(invKbb)*Dbloc';
        [jLloc,iLloc] = meshgrid(dofId(i,3),dofId(i,3));
        nl = numel(Lloc);
        iLvec(cl+1:cl+nl) = iLloc(:);   jLvec(cl+1:cl+nl) = jLloc(:);
        Lvec(cl+1:cl+nl) = Lloc(:);

        % update counters
        cd = cd+nd;
        ck = ck+nk;
        cl = cl+nl;
        cb = cb+nb;
         
        % track element not fully projected
        if ~all(id)
          error('%i GP not sorted for slave elem numb %i \n',sum(id),is);
        end
      end

      % cut master index vectors for sparse matrix assembly
      iMvec = iMvec(1:cm); jMvec = jMvec(1:cm); Mvec = Mvec(1:cm);

      % assemble mortar matrices in sparse format
      obj.Jmaster{1} = sparse(iMvec,jMvec,-Mvec,nDofMult,nDofMaster);
      obj.Jslave{1} = sparse(iDvec,jDvec,Dvec,nDofMult,nDofSlave)+...
        sparse(iBvec,jBvec,Bvec,nDofMult,nDofSlave);
      obj.Jmult{1} = sparse(iLvec,jLvec,Lvec,nDofMult,nDofMult); 
%       % apply condensation to slave inner block
      Jcondensation = sparse(iKvec,jKvec,Kvec,nDofSlave,nDofSlave);
      poro.J = poro.J + Jcondensation;
    end

    function computeMat(obj,dt)
      computeMortarMatrices(obj,dt);
    end


    function out = isMatrixComputed(obj)
      out = all(cellfun(@(x) ~isempty(x), ...
        [obj.Jmaster(:); obj.Jslave(:); obj.Jmult(:)]));
    end

    function updateState(obj,du)
      % enhance strain with bubble contribution

      for i = 1:obj.nFld
        n = numel(obj.multipliers(i).curr);
        obj.multipliers(i).curr = obj.multipliers(i).curr + du(1:n);
        du = du(n+1:end);
        % update also strain considering the bubble finite element
        % this allow to compute the rhs using the classical Bu^T sigma
      end
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

