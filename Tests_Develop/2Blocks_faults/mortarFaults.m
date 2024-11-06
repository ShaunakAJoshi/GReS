classdef mortarFaults < handle
   properties
      simParameters
      t = 0
      nDom = 2
      nInterf = 1
      tStep = 0
      iter
      dt
      state
      models
      meshGlue
      Dg   % mortar slave matrix in global coords
      Mg   % mortar cross matrix in global coords
      Dn   % mortar slave matrix for normal gaps
      Mn   % mortar cross matrix for normal gaps
      Dt   % mortar slave matrix for tangential gaps
      Mt   % mortar cross matrix for tangential gaps
      nNslave
      nNmaster
      indN % index of basis functions in displacement basis matrix
      nn_s
      nn_m
   end

   methods
      function obj = mortarFaults(simParam,meshGlue)
         obj.setParams(simParam,meshGlue);
         obj.simulationLoop();
      end

      function simulationLoop(obj)
         obj.dt = obj.simParameters.dtIni;
         delta_t = obj.dt; % dynamic time step
         % compute mechanical matrices
         for i = 1:2
            computeLinearMatrices(obj.meshGlue.model(i).Discretizer,obj.state(i).curr,obj.state(i).prev,obj.dt)
         end
         flConv = true; %convergence flag
         % Loop over time
         while obj.t < obj.simParameters.tMax
            % Update the simulation time and time step ID
            absTol = obj.simParameters.absTol;
            obj.tStep = obj.tStep + 1;
            %new time update to fit the outTime list
            [obj.t, delta_t] = obj.updateTime(flConv, delta_t);
            for i = 1:obj.nDom
               obj.meshGlue.model(i).Discretizer.resetJacobianAndRhs();
               % Apply the Dirichlet condition value to the solution vector
               if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                  applyDirVal(obj.meshGlue.model(i).ModelType,obj.meshGlue.model(i).BoundaryConditions,...
                     obj.t, obj.state(i).curr);
               end
               %
            end

            % Conctact mechanics algorithms
            getContactMatrices(obj);




            % compute Rhs norm
            rhsNorm = norm(rhs,2);

            if obj.simParameters.verbosity > 0
               fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
               fprintf('-----------------------------------------------------------\n');
            end

            if obj.simParameters.verbosity > 1
               fprintf('Iter     ||rhs||\n');
            end

            tolWeigh = obj.simParameters.relTol*rhsNorm;
            obj.iter = 0;
            %
            if obj.simParameters.verbosity > 1
               fprintf('0     %e\n',rhsNorm);
            end

            while ((rhsNorm > tolWeigh) && (obj.iter < obj.simParameters.itMaxNR) ...
                  && (rhsNorm > absTol)) || obj.iter == 0
               obj.iter = obj.iter + 1;
               %
               du = J\-rhs;
               % update solution vector for each model
               obj.updateStateMD(du);
               for i = 1:obj.nDom
                  obj.meshGlue.model(i).Discretizer.resetJacobianAndRhs();
                  % Update tmpState
                  computeNLMatricesAndRhs(obj.meshGlue.model(i).Discretizer,...
                     obj.state(i).curr,obj.state(i).prev,obj.dt);
                  % compute block Jacobian and block Rhs
                  obj.meshGlue.model(i).Discretizer.computeBlockJacobianAndRhs(delta_t);
               end


               % Get unique multidomain solution system
               [J,rhs] = obj.meshGlue.getMDlinSyst();

               % Apply BCs to global linear system
               for i = 1:obj.nDom
                  if ~isempty(obj.meshGlue.model(i).BoundaryConditions)
                     [J,rhs] = applyBCAndForces_MD(i, obj.meshGlue, obj.t, obj.state(i).curr, J, rhs);
                  end
               end

               % compute Rhs norm
               rhsNorm = norm(rhs,2);
               if obj.simParameters.verbosity > 1
                  fprintf('%d     %e\n',obj.iter,rhsNorm);
               end
            end
            %
            % Check for convergence
            flConv = (rhsNorm < tolWeigh || rhsNorm < absTol);
            if flConv % Convergence
               for i = 1:obj.nDom
                  obj.state(i).curr.t = obj.t;
                  % Print the solution, if needed
                  if isPoromechanics(obj.meshGlue.model(i).ModelType)
                     obj.state(i).curr.advanceState();
                  end
                  if isVariabSatFlow(obj.meshGlue.model(i).ModelType)
                     obj.state(i).curr.updateSaturation()
                  end
               end
               if obj.t > obj.simParameters.tMax   % For Steady State
                  for i = 1:obj.nDom
                     printState(obj.meshGlue.model(i).OutState,obj.state(i).curr);
                  end
               else
                  for i = 1:obj.nDom
                     printState(obj.meshGlue.model(i).OutState,obj.state(i).prev,obj.state(i).curr);
                  end
               end
            end

            %
            % Manage next time step
            delta_t = manageNextTimeStep(obj,delta_t,flConv);
         end
         %
      end

      function setParams(obj,simParam,mG)
         obj.simParameters = simParam;
         obj.models = mG.model;
         obj.meshGlue = mG;
         obj.nNmaster = mG.interfaces.mortar.nNmaster;
         obj.nNslave = mG.interfaces.mortar.nNslave;
         assert(numel(obj.models)==2,'Too many input domains');
         % state structure easier to access
         obj.state = repmat(struct('prev',{},'curr',{}),2,1);
         obj.state(1).prev = obj.meshGlue.model(1).State;
         obj.state(1).curr = copy(obj.state(1).prev);
         obj.state(2).prev = obj.meshGlue.model(2).State;
         obj.state(2).curr = copy(obj.state(2).prev);
         % pattern of basis functions transform
         obj.indN = repmat([1 5 9],1,obj.nNslave)+9*repelem(0:obj.nNslave-1,1,3);
         obj.nn_s = length(mG.interfaces.mortar.nodesSlave);
         obj.nn_m = length(mG.interfaces.mortar.nodesMaster);
      end

      function [dt] = manageNextTimeStep(obj,dt,flConv)
         if ~flConv   % Perform backstep
            for i = 1:obj.nDom
               transferState(obj.state(i).curr,obj.state(i).prev);
            end
            obj.t = obj.t - obj.dt;
            obj.tStep = obj.tStep - 1;
            dt = dt/obj.simParameters.divFac;
            obj.dt = obj.dt/obj.simParameters.divFac;  % Time increment chop
            if min(dt,obj.dt) < obj.simParameters.dtMin
               if obj.simParameters.goOnBackstep == 1
                  flConv = 1;
               elseif obj.simParameters.goOnBackstep == 0
                  error('Minimum time step reached')
               end
            elseif obj.simParameters.verbosity > 0
               fprintf('\n %s \n','BACKSTEP');
            end
         end
         if flConv % Go on if converged
            obj.dt = min([obj.dt*obj.simParameters.multFac,obj.simParameters.dtMax]);
            obj.dt = max([obj.dt obj.simParameters.dtMin]);
            for i = 1:obj.nDom
               transferState(obj.state(i).curr,obj.state(i).prev);
            end
            %
            if ((obj.t + obj.dt) > obj.simParameters.tMax)
               obj.dt = obj.simParameters.tMax - obj.t;
            end
         end
      end


      function [t, dt] = updateTime(obj,conv,dt)
         tMax = obj.simParameters.tMax;
         for i = 1:obj.nDom
            if obj.meshGlue.model(i).OutState.modTime
               tmp = find(obj.t<obj.meshGlue.model(i).outState.timeList(),1,'first');
               if ~conv
                  t = min([obj.t + obj.dt, obj.t + dt, obj.meshGlue.model(i).OutState.timeList(tmp)]);
               else
                  t = min([obj.t + obj.dt, obj.meshGlue.model(i).OutState.timeList(tmp)]);
               end
            else
               t = obj.t + obj.dt;
            end
            if t > tMax
               t = tMax;
            end
         end
         dt = t - obj.t;
      end


      function getContactMatrices(obj)
         % get and normalize nodal normals
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         n = obj.meshGlue.interfaces(1).nodeNormal./area_nod;
         % compute mortar matrices within one simulation loop
         % use radial basis functions to evaluate matrix M
         mortar = obj.meshGlue.interfaces(1).mortar; % only one interface defined
         nGP = 4; % numb. of gauss points for element based integration
         nInt = 4; % numb. of interpolation points for RBF intrpolation
         tol = 1e-3;
         type = 'gauss';
         mult_type = 'dual';
         c_ns = 0;  % counter for GP not projected
         %Mdetect = zeros(mortar.nElMaster,mortar.nElSlave);
         % set Gauss class
         gM = Gauss(mortar.masterCellType,3,2); % gauss class for Master element interpolation
         gS = Gauss(mortar.slaveCellType,nGP,2); % gauss class for slave integration
         elemMaster = getElem(mortar,gM,'master');
         elemSlave = getElem(mortar,gS,'slave');
         [imVec,jmVec,Mgvec,Mtvec,Mnvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         [isVec,jsVec,Dgvec,Dtvec,Dnvec] = deal(zeros(nnz(mortar.elemConnectivity)*mortar.nNmaster^2,1));
         % Perform interpolation on the master side (computing weights
         % and interpolation coordinates)
         [wFMat,w1Mat,ptsIntMat] = mortar.getWeights('master',nInt,elemMaster,type);
         %[wFMatS,w1MatS,ptsIntMatS] = getWeights(obj,'slave',nInt,elemSlave,type);
         % Interpolation for support detection
         %[wFSupp,w1Supp] = getSuppWeight(obj);
         % Loop trough slave elements
         cs = 0; % slave matrix entry counter
         cm = 0; % master matrix entry counter
         for j = 1:mortar.nElSlave
            %Compute Slave quantities
            dJWeighed = elemSlave.getDerBasisFAndDet(j,3); % Weighted Jacobian
            %get Gauss Points position in the real space
            ptsGauss = getGPointsLocation(elemSlave,j);
            idSlave = mortar.intSlave.surfaces(j,:);
            n_el = (n(idSlave,:))';
            n_el = n_el(:);
            R = getRotationMatrix(obj,n_el);
            A = diag(repelem(area_nod(idSlave),3));
            NSlave = getBasisFinGPoints(elemSlave); % Slave basis functions
            switch mult_type
               case 'standard'
                  NSlaveMult = NSlave; % Slave basis functions
               case 'dual'
                  NSlaveMult = mortar.computeDualBasisF(NSlave,dJWeighed);
            end
            master_elems = find(mortar.elemConnectivity(:,j));
            for jm = master_elems'
               idMaster = mortar.intMaster.surfaces(jm,:);
               ptsInt = ptsIntMat(:,repNum(3,jm));
               [fiNM,id1] = mortar.computeRBFfiNM(ptsInt,ptsGauss,type);
               switch mortar.degree
                  case 1
                     NMaster = (fiNM*wFMat(:,repNum(mortar.nNmaster,jm)))./(fiNM*w1Mat(:,jm));
                     Nsupp = NMaster(:,[1 2 3]);
                  case 2
                     Ntmp = (fiNM*wFMat(:,repNum(mortar.nNmaster+2,jm)))./(fiNM*w1Mat(:,jm));
                     NMaster = Ntmp(:,1:mortar.nNmaster);
                     Nsupp = Ntmp(:,[end-1 end]);
               end
               % automatically detect supports computing interpolant
               id = all([Nsupp >= 0-tol id1],2);
               if any(id)
                  % element-based integration
                  % prepare provisional 3D matrices
                  Nm = obj.dispSP(NMaster(id,:));
                  Ns = obj.dispSP(NSlave(id,:));
                  Nmult = obj.dispSP(NSlaveMult(id,:));
                  % global mortar matrices
                  Dgtmp = pagemtimes(Nmult,'transpose',Ns,'none');
                  Dgtmp = Dgtmp.*reshape(dJWeighed(id),1,1,[]);
                  Dgloc = R'*sum(Dgtmp,3);
                  Mgtmp = pagemtimes(Nmult,'transpose',Nm,'none');
                  Mgtmp = Mgtmp.*reshape(dJWeighed(id),1,1,[]);
                  Mgloc = R'*sum(Mgtmp,3);
                  % normal mortar matrices
                  Nn = pagemtimes(Ns,n_el);
                  Dntmp = pagemtimes(Nmult(1,:,:),'transpose',pagemtimes(Nn,'transpose',Ns,'none'),'none');
                  Dntmp = Dntmp.*reshape(dJWeighed(id),1,1,[]);
                  Dnloc = sum(Dntmp,3);
                  Mntmp = pagemtimes(Nmult(1,:,:),'transpose',pagemtimes(Nn,'transpose',Nm,'none'),'none');
                  Mntmp = Mntmp.*reshape(dJWeighed(id),1,1,[]);
                  Mnloc = sum(Mntmp,3);
                  % tangential mortar matrices
                  Nt = repmat(eye(3),[1,1,[]]) - pagemtimes(Nn,'none',Nn,'transpose');
                  Dttmp = pagemtimes(Nmult,'transpose',pagemtimes(Nt,Ns),'none');
                  Dttmp = Dttmp.*reshape(dJWeighed(id),1,1,[]);
                  Dtloc = A*sum(Dttmp,3);
                  Mttmp = pagemtimes(Nmult,'transpose',pagemtimes(Nt,Nm),'none');
                  Mttmp = Mttmp.*reshape(dJWeighed(id),1,1,[]);
                  Mtloc = A*sum(Mttmp,3);
                  dof_master = obj.get_dof(idMaster);
                  dof_slave = obj.get_dof(idSlave);
                  [jjM,iiM] = meshgrid(dof_master,dof_slave);
                  [jjS,iiS] = meshgrid(dof_slave,dof_slave);
                  nm = numel(Mgloc);
                  ns = numel(Dgloc);
                  imVec(cm+1:cm+nm) = iiM(:); jmVec(cm+1:cm+nm) = jjM(:);
                  isVec(cs+1:cs+ns) = iiS(:); jsVec(cs+1:cs+ns) = jjS(:);
                  Mgvec(cm+1:cm+nm) = Mgloc(:);
                  Dgvec(cs+1:cs+ns) = Dgloc(:);
                  Mnvec(cm+1:cm+nm) = Mnloc(:);
                  Dnvec(cs+1:cs+ns) = Dnloc(:);
                  Mtvec(cm+1:cm+nm) = Mtloc(:);
                  Dtvec(cs+1:cs+ns) = Dtloc(:);
                  % sort out Points already projected
                  dJWeighed = dJWeighed(~id);
                  ptsGauss = ptsGauss(~id,:);
                  NSlave = NSlave(~id,:);
                  NSlaveMult = NSlaveMult(~id,:);
                  cs = cs+ns;
                  cm = cm+nm;
               end
            end
            if ~all(id)
               fprintf('GP not sorted for slave elem %i \n',j);
               c_ns = c_ns + 1;
            end
         end
         imVec = imVec(1:cm); jmVec = jmVec(1:cm);
         isVec = isVec(1:cs); jsVec = jsVec(1:cs);
         Mgvec = Mgvec(1:cm); Mnvec = Mnvec(1:cm); Mtvec = Mtvec(1:cm);
         Dgvec = Dgvec(1:cm); Dnvec = Dnvec(1:cm); Dtvec = Dtvec(1:cm);
         obj.Mg = sparse(imVec,jmVec,Mgvec,3*obj.nn_s,3*obj.nn_m);
         obj.Mn = sparse(imVec,jmVec,Mnvec,3*obj.nn_s,3*obj.nn_m);
         obj.Mt = sparse(imVec,jmVec,Mtvec,3*obj.nn_s,3*obj.nn_m);
         obj.Dg = sparse(isVec,jsVec,Dgvec,3*obj.nn_s,3*obj.nn_s);
         obj.Dn = sparse(isVec,jsVec,Dnvec,3*obj.nn_s,3*obj.nn_s);
         obj.Dt = sparse(isVec,jsVec,Dtvec,3*obj.nn_s,3*obj.nn_s);
      end

      function R = getRotationMatrix(obj,n_el)
         % n_el: local nodal normal array
         R = zeros(3*obj.nNslave,3*obj.nNslave);
         for i = 0:obj.nNslave-1
            n = n_el(3*i+1:3*i+3);
            if all(abs(n) ~= [1;0;0])
               u = [1; 0; 0];
            else
               u = [0; 1; 0];
            end
            m1 = cross(n,u)./norm(cross(n,u),2);
            m2 = cross(n,m1);
            R(3*i+1:3*i+3,3*i+1:3*i+3) = [n m1 m2];
         end
      end

      function Nout = dispSP(obj,Nin)
         % reshape basis functions to obtain displacement shape function
         % input: nG x nN matrix
         % output: 3 x nN x nG
         s2 = 3*obj.nNslave;
         s3 = size(Nin,1);
         N = zeros(3*s2,s3);
         Nin = repelem(Nin',3,1);
         N(obj.indN,:) = Nin;
         Nout = reshape(N,[3,s2,s3]); % reshaped 3D matrix
      end
   end

   methods (Static)
      function dof = get_dof(slaveNodes)
         dof = repelem(3*slaveNodes,3);
         dof = dof + repmat(-2:0,1,numel(slaveNodes));
      end
   end
end

