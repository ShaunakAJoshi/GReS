classdef mortarFaults < handle
   properties
      nDom = 2
      nInterf = 1
      meshGlue    % instance of the meshGlue class (contains mortar algorithms)
      % MORTAR MATRICES
      Dg          % mortar slave matrix in global coords
      Mg          % mortar cross matrix in global coords
      Dn          % mortar slave matrix for normal gaps
      Mn          % mortar cross matrix for normal gaps
      Dt          % mortar slave matrix for tangential gaps
      Mt          % mortar cross matrix for tangential gaps
      L           % lagrange multiplier mass matrix (local coordinates)
      nNslave     % number of nodes per slave elements
      nNmaster    % number
      indN        % index of basis functions in displacement basis matrix
      nS          % number of interface master nodes
      nM          % number of interface slave nodes
      coes        % coesion
      phi         % friction angle
      tagMaster   % id of master domain
      tagSlave     % id of slave domain
      idMaster
      idSlave
      totNodMaster % total number of master domain nodes
      totNodSlave  % total number of slave domain nodes
      nDoF
      activeSet   % group nodes based on active set state
      gap         % g
      g_T         % \Delta_n g_T
      g_N         % n \cdot g
      E           % mortar operator
      R           % global rotation matrix
      tolGap = 1e-5; % tolerance on normal and tangential gaps
      tolNormal = 1e-3 % tolerance on normal traction
      tolTang = 1e-2 % relative tolerance w.r.t tau_lim
   end

   methods
      function obj = mortarFaults(meshGlue,coes,phi)
         obj.coes = coes;
         obj.phi = phi;
         obj.setParams(meshGlue);
      end

      function setParams(obj,mG)
         obj.meshGlue = mG;
         obj.nNmaster = mG.interfaces.mortar.nNmaster;
         obj.nNslave = mG.interfaces.mortar.nNslave;
         assert(numel(mG.model)==2,'Too many input domains');
         % pattern of basis functions transform
         obj.indN = repmat([1 5 9],1,obj.nNslave)+9*repelem(0:obj.nNslave-1,1,3);
         % save model metrics
         obj.tagMaster = obj.meshGlue.interfaces(1).Master;
         obj.tagSlave = obj.meshGlue.interfaces(1).Slave;
         obj.totNodMaster = obj.meshGlue.model(obj.tagMaster).Grid.topology.nNodes;
         obj.totNodSlave = obj.meshGlue.model(obj.tagSlave).Grid.topology.nNodes;
         obj.nS = length(mG.interfaces.mortar.nodesSlave);
         obj.nM = length(mG.interfaces.mortar.nodesMaster);
         obj.idMaster = obj.meshGlue.interfaces(1).masterSet;
         obj.idSlave = obj.meshGlue.interfaces(1).slaveSet;
         % compute global rotation matrix
         obj.R = getGlobalRotationMatrix(obj);
         % initialize stress (not flexible at all)
         % build dof map with initial active set status
         %
         obj.nDoF = 3*obj.totNodMaster+3*obj.totNodSlave+3*obj.nS;
         obj.E = obj.meshGlue.interfaces.mortar.getMortarOperator(3); % mortar operator
      end


      function computeContactMatrices(obj)
         % get and normalize nodal normals
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2)); % nodal area
         n = obj.meshGlue.interfaces(1).nodeNormal./area_nod;
         % compute mortar matrices within one simulation loop use radial
         % basis functions to evaluate matrix M
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
            nSlave = mortar.intSlave.surfaces(j,:);
            n_el = (n(nSlave,:))';
            n_el = n_el(:); % local nodal normal vector
            Rloc = -getRotationMatrix(obj,n_el);
            %Rloc = eye(3*obj.nNslave);
            %A = diag(repelem(area_nod(nSlave),3));
            NSlave = getBasisFinGPoints(elemSlave); % Slave basis functions
            switch mult_type
               case 'standard'
                  NSlaveMult = NSlave; % Slave basis functions
               case 'dual'
                  NSlaveMult = mortar.computeDualBasisF(NSlave,dJWeighed);
            end
            master_elems = find(mortar.elemConnectivity(:,j));
            for jm = master_elems'
               nMaster = mortar.intMaster.surfaces(jm,:);
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
                  % global mortar matrices (global frame)
                  Dgtmp = pagemtimes(Nmult,'transpose',Ns,'none');
                  Dgtmp = Dgtmp.*reshape(dJWeighed(id),1,1,[]);
                  Dgloc = Rloc'*sum(Dgtmp,3);
                  Mgtmp = pagemtimes(Nmult,'transpose',Nm,'none');
                  Mgtmp = Mgtmp.*reshape(dJWeighed(id),1,1,[]);
                  Mgloc = Rloc'*sum(Mgtmp,3);
                  % normal mortar matrices (global frame)
                  Nn = pagemtimes(Ns,n_el);
                  Dntmp = pagemtimes(Nmult(1,:,:),'transpose',pagemtimes(Nn,'transpose',Ns,'none'),'none');
                  Dntmp = Dntmp.*reshape(dJWeighed(id),1,1,[]);
                  Dnloc = Rloc'*sum(Dntmp,3);
                  Mntmp = pagemtimes(Nmult(1,:,:),'transpose',pagemtimes(Nn,'transpose',Nm,'none'),'none');
                  Mntmp = Mntmp.*reshape(dJWeighed(id),1,1,[]);
                  Mnloc = Rloc'*sum(Mntmp,3);
                  % tangential mortar matrices (global frame)
                  Nt = repmat(eye(3),[1,1,[]]) - pagemtimes(Nn,'none',Nn,'transpose');
                  Dttmp = pagemtimes(Nmult,'transpose',pagemtimes(Nt,Ns),'none');
                  Dttmp = Dttmp.*reshape(dJWeighed(id),1,1,[]);
                  Dtloc = Rloc'*sum(Dttmp,3);
                  Mttmp = pagemtimes(Nmult,'transpose',pagemtimes(Nt,Nm),'none');
                  Mttmp = Mttmp.*reshape(dJWeighed(id),1,1,[]);
                  Mtloc = Rloc'*sum(Mttmp,3);
                  % lagrange multipliers mass matrix L (local frame!)
                  Ltmp = pagemtimes(Nmult,'transpose',Nmult,'none');
                  Ltmp = Ltmp.*reshape(dJWeighed(id),1,1,[]);
                  Lloc = sum(Ltmp,3);
                  % local assembly
                  dof_master = get_dof(nMaster);
                  dof_slave = get_dof(nSlave);
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
                  Lvec(cs+1:cs+ns)  = Lloc(:);
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
         Lvec = Lvec(1:cs);
         obj.Mg = sparse(imVec,jmVec,Mgvec,3*obj.nS,3*obj.nM);
         obj.Mn = sparse(imVec,jmVec,Mnvec,3*obj.nS,3*obj.nM);
         obj.Mt = sparse(imVec,jmVec,Mtvec,3*obj.nS,3*obj.nM);
         obj.Dg = sparse(isVec,jsVec,Dgvec,3*obj.nS,3*obj.nS);
         obj.Dn = sparse(isVec,jsVec,Dnvec,3*obj.nS,3*obj.nS);
         obj.Dt = sparse(isVec,jsVec,Dtvec,3*obj.nS,3*obj.nS);
         obj.L = sparse(isVec,jsVec,Lvec,3*obj.nS,3*obj.nS);
         % eliminate small numerical quantities in obj.Dg
         obj.Dg(abs(obj.Dg)<1e-12) = 0;
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

      function R = getNodeRotationMatrix(obj,k)
         % n_el: local nodal normal array
         n = obj.meshGlue.interfaces.nodeNormal(k,:);
         n = (n/norm(n,2))';
         if all(abs(n) ~= [1;0;0])
            u = [1; 0; 0];
         else
            u = [0; 1; 0];
         end
         m1 = cross(n,u)./norm(cross(n,u),2);
         m2 = cross(n,m1);
         R = [n m1 m2];
      end
      

      function R = getGlobalRotationMatrix(obj)
         area_nod = sqrt(sum(obj.meshGlue.interfaces(1).nodeNormal.^2,2));
         n_el = obj.meshGlue.interfaces(1).nodeNormal./area_nod; % normalized nodal area
         % n_el: local nodal normal array
         nEntry = 9;
         Rvec = zeros(nEntry*obj.nS,1);
         iivec = zeros(nEntry*obj.nS,1);
         jjvec = zeros(nEntry*obj.nS,1);
         l1 = 0;
         for i = 0:obj.nS-1
            n = n_el(i+1,:);
            if all(abs(n) ~= [1;0;0])
               u = [1; 0; 0];
            else
               u = [0; 1; 0];
            end
            m1 = cross(n,u)./norm(cross(n,u),2);
            m2 = cross(n,m1);
            [iiLoc,jjLoc] = meshgrid(3*i+1:3*i+3,3*i+1:3*i+3);
            Rloc = [n' m1' m2'];
            iivec(l1+1:l1+nEntry) = iiLoc(:);
            jjvec(l1+1:l1+nEntry) = jjLoc(:);
            Rvec(l1+1:l1+nEntry) = Rloc(:);
            l1 = l1+nEntry;
         end
         R = sparse(iivec,jjvec,Rvec,3*obj.nS,3*obj.nS);
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

      function dofs = getContactDofs(obj,tag,varargin)
         if ~isempty(varargin)
            dofIn = varargin{1};
            dofIn = reshape(dofIn,1,[]);
         else
            dofIn = [];
         end
         % dofIn must be a row vector
         % map nodal entries to degrees of freedom in the global linear
         % system. 
         % tag: 'master','slave','lag'
         switch tag
            case 'master'
               if isempty(dofIn)
                  dofIn = 1:obj.totNodMaster;
               end
               dofs = get_dof(dofIn);
            case 'slave'
               if isempty(dofIn)
                  dofIn = 1:obj.totNodSlave;
               end
               dofs = 3*obj.totNodMaster+get_dof(dofIn);
            case 'lag'
               if isempty(dofIn)
                  dofs = [];
                  return
               end
               dofs = 3*(obj.totNodMaster+obj.totNodSlave)+get_dof(dofIn);
            otherwise
               error(['Invalide tag string for getContactDofs method \n:' ...
                  'valid inputs are: master, slave, lag']);
         end
      end

      function tLim = computeLimitTraction(obj,activeSet)
         % TO DO: add check on magnitud gap
         % get limit traction array in local coordinates
         tLim = zeros(3*numel(activeSet.new.slip),1);
         for i = activeSet.new.slip'
            dofSlip = get_dof(i);
            gapNorm = norm(obj.g_T(i),2);
            if gapNorm < obj.tolGap
               tL = zeros(3,1);
            else
               tL = repelem(tau_max(obj,i),3,1).*obj.g_T(dofSlip)/gapNorm; % traction in global coords
               tL = obj.getNodeRotationMatrix(i)'*tL; % traction in local coords
               tL(1) = 0; % set normal component to 0 for rhs computation
            end
            tLim(3*i-2:3*i) = tL;
         end
      end

      function tau = tau_max(obj,i,multipliers)
         % use the normal component of contact traction
         tau = obj.coes - tan(deg2rad(obj.phi))*multipliers(3*i-1);
      end

      function computeNodalGap(obj,state,dofMap)
         % compute gap in global coordinates
         % compute time difference of the tangential component of nodal gap
         % \Delta_n g_T (in global coordinates)
         % Use mortar operator E to map master nodes to slave side
         usCurr = state(obj.tagSlave).dispCurr(dofMap.nodSlave);
         umCurr = state(obj.tagMaster).dispCurr(dofMap.nodMaster);
         g_old = obj.gap;
         obj.gap = usCurr - obj.E*umCurr;
         obj.g_T = obj.gap - g_old; % global nodal gap
         area_nod = sqrt(sum(obj.meshGlue.interfaces.nodeNormal.^2,2)); % nodal area
         n = obj.meshGlue.interfaces.nodeNormal./area_nod; % unit nodal normal
         % compute tangential component
         obj.g_N = zeros(obj.nS,1);
         l1 = 0;
         for i = 1:obj.nS
            % get node normal
            n_i = n(i,:);
            T = eye(3) - n_i'*n_i; % tangential projection matrix
            obj.g_T(l1+1:l1+3) = T*obj.g_T(l1+1:l1+3);
            obj.g_N(i) = n_i*obj.gap(l1+1:l1+3);
            l1 = l1+3;
         end
      end

      function Tmat = computeDtDgt(obj,activeSet)
         % compute derivative of tangential traction w.r.t tangential gap
         % quantities are in global coordinates
         iVec = zeros(9*numel(activeSet.curr.slip),1);
         jVec = zeros(9*numel(activeSet.curr.slip),1);
         Tvec = zeros(9*numel(activeSet.curr.slip),1);
         l = 0;
         for i = activeSet.curr.slip'
            dof = get_dof(i);
            tLim = tau_max(obj,i);
            DgT = obj.g_T(dof);
            normgT = norm(DgT);
            if normgT > obj.tolGap
               Tloc = tLim*(normgT^2*eye(3)-Dgt*Dgt')/normgT^3; % 3x3 local mat in global coords
               [ii,jj] = meshgrid(dof,dof);
               iVec(l+1:l+9) = ii(:);
               jVec(l+1:l+9) = jj(:);
               Tvec(l+1:l+9) = Tloc(:);
               l = l+9;
            end
         end
         Tmat = sparse(iVec,jVec,Tvec,3*obj.nS,3*obj.nS);
      end

      function Nmat = computeDtDtn(obj,activeSet)
         % compute derivative of tangential traction w.r.t tangential gap
         % quantities are in local coordinates
         iVec = zeros(2*numel(activeSet.curr.slip),1);
         jVec = zeros(2*numel(activeSet.curr.slip),1);
         Nvec = zeros(2*numel(activeSet.curr.slip),1);
         l = 0;
         for i = activeSet.curr.slip'
            dof = get_dof(i);
            tLim = tau_max(obj,i);
            DgT = obj.getNodeRotationMatrix(i)'*obj.g_T(dof); % local tangential gap
            DgT = DgT([2 3]); % if evrything is ok, the first component is actually 0
            normgT = norm(DgT);
            if normgT > obj.tolGap
               Nloc = tLim*gT/normgT; % 3x3 local mat in global coords
               iVec(l+1:l+2) = [dof(2);dof(3)];
               jVec(l+1:l+2) = [dof(1);dof(1)];
               Nvec(l+1:l+2) = Nloc(:);
               l = l+2;
            end
         end
         Nmat = sparse(iVec,jVec,Nvec,3*obj.nS,3*obj.nS);
      end
   end
end
