classdef MultidomainFCSolver < handle
  % Class for solving non linear problem involving multiple non conforming
  % domains using the mortar method

  properties (Access = private)
    %
    simParameters
    nDom
    nInterf
    %
    t = 0
    tStep = 0
    iter
    dt
    statek
    stateTmp
    systSize
    nfldInt
    nfldDom
  end

  properties (Access = public)
    state
    domains
    interfaces
  end

  methods (Access = public)
    function obj = MultidomainFCSolver(simParam,models,interfaces)
      obj.setNonLinearSolver(simParam,models,interfaces);
    end

    function NonLinearLoop(obj)

      % Initialize the time step increment
      obj.dt = obj.simParameters.dtIni;
      delta_t = obj.dt; % dynamic time step

      %
      flConv = true; %convergence flag
      %

      % Loop over time
      while obj.t < obj.simParameters.tMax
        % Update the simulation time and time step ID
        absTol = obj.simParameters.absTol;
        obj.tStep = obj.tStep + 1;
        %new time update to fit the outTime list

        [obj.t, delta_t] = obj.updateTime(flConv, delta_t);

        if obj.simParameters.verbosity > 0
          fprintf('\nTSTEP %d   ---  TIME %f  --- DT = %e\n',obj.tStep,obj.t,delta_t);
          fprintf('-----------------------------------------------------------\n');
        end
        if obj.simParameters.verbosity > 1
          fprintf('Iter     ||rhs||\n');
        end

        for i = 1:obj.nDom

          discretizer = obj.domains(i).Discretizer;
          bc = obj.domains(i).BoundaryConditions;

          % Check if boundary conditions are defined for the i-th domain
          if ~isempty(obj.domains(i).BoundaryConditions)
            % Apply Dirichlet boundary values to i-th domain
            obj.state(i).curr = applyDirVal(...
              discretizer, bc, obj.t, obj.state(i).curr);
          end

          % Compute matrices and residuals for individual models of
          % the i-th domain
          obj.state(i).curr = computeMatricesAndRhs(...
            discretizer, obj.state(i).curr, obj.state(i).prev, obj.dt);

          % Compute domain coupling matrices and rhs
          for j = discretizer.interfaceList
            computeMat(obj.interfaces{j},i);
            computeRhs(obj.interfaces{j},i,obj.state(i).curr);
          end

          % Apply BCs to the blocks of the linear system
          applyBC(discretizer, bc, obj.t, obj.state(i).curr);

          % Apply BC to coupling matrices
          for j = discretizer.interfaceList
            applyBC(obj.interfaces{j},i,bc,obj.t, obj.state(i).curr);
          end

        end

        % Assemble multidomain solution system
        rhs = assembleRhs(obj);

        % compute Rhs norm
        rhsNorm = norm(cell2mat(rhs),2);

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
          J = assembleJacobian(obj);

          du = FCSolver.solve(J,rhs);
          % update solution vector for each model
          obj.updateStateMD(du);
          for i = 1:obj.nDom
            obj.domains(i).Discretizer.resetJacobianAndRhs();
            % Update tmpState
            computeNLMatricesAndRhs(obj.domains(i).Discretizer,...
              obj.state(i).curr,obj.state(i).prev,obj.dt);
            % compute block Jacobian and block Rhs
            obj.domains(i).Discretizer.computeBlockJacobianAndRhs(delta_t);
          end


          % Get unique multidomain solution system
          [J,rhs] = obj.meshGlue.getMDlinSyst();

          % Apply BCs to global linear system
          for i = 1:obj.nDom
            if ~isempty(obj.domains(i).BoundaryConditions)
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
            if isPoromechanics(obj.domains(i).ModelType)
              obj.state(i).curr.advanceState();
            end
            if isVariabSatFlow(obj.domains(i).ModelType)
              obj.state(i).curr.updateSaturation()
            end
          end
          if obj.t > obj.simParameters.tMax   % For Steady State
            for i = 1:obj.nDom
              printState(obj.domains(i).OutState,obj.state(i).curr);
            end
          else
            for i = 1:obj.nDom
              printState(obj.domains(i).OutState,obj.state(i).prev,obj.state(i).curr);
            end
          end
        end

        %
        % Manage next time step
        delta_t = manageNextTimeStep(obj,delta_t,flConv);
      end
      %
    end
  end

  methods (Access = private)
    function setNonLinearSolver(obj,simParam,dom,interf)
      obj.simParameters = simParam;
      obj.domains = dom;
      obj.nDom = numel(dom);
      obj.interfaces = interf;
      obj.nInterf = numel(interf);
      obj.state = repmat(struct('prev',{},'curr',{}),obj.nDom,1);
      % initialize a state structure for each domains
      for i = 1:obj.nDom
        obj.state(i).prev = obj.domains(i).State;
        obj.state(i).curr =  obj.domains(i).State;
      end
      getSystemSize(obj);
      getNumField(obj);
    end

    function [t, dt] = updateTime(obj,conv,dt)
      t = obj.simParameters.tMax;
      told = t;
      for i = 1:obj.nDom
        if obj.domains(i).OutState.modTime
          tmp = find(obj.t<obj.domains(i).outState.timeList(),1,'first');
          if ~conv
            t = min([obj.t + obj.dt, obj.t + dt, obj.domains(i).OutState.timeList(tmp)]);
          else
            t = min([obj.t + obj.dt, obj.domains(i).OutState.timeList(tmp)]);
          end
        else
          t = obj.t + obj.dt;
        end
        if t > told
          t = told;
        end
      end
      dt = t - obj.t;
    end

    function out = computeRhsNorm(obj)
      %Return maximum norm of the entire domain
      rhsNorm = zeros(obj.nDom,1);
      for i = 1:obj.nDom
        nRhs = length(obj.domains(i).DoFManager.subList);
        rhsNorm_loc = zeros(nRhs,1);
        for j = 1:nRhs
          rhsNorm_loc(j) = norm(obj.domains(i).Discretizer.rhs{j}, obj.simParameters.pNorm);
        end
        rhsNorm(i) = sqrt(sum(rhsNorm_loc.^2));
      end
      out = norm(rhsNorm);
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
        tmpVec = obj.simParameters.multFac;
        for i =1:obj.nDom
          if isFlow(obj.domains(i).ModelType)
            dpMax = max(abs(obj.state(i).curr.pressure - obj.state(i).prev.pressure));
            tmpVec = [tmpVec, (1+obj.simParameters.relaxFac)* ...
              obj.simParameters.pTarget/(dpMax + obj.simParameters.relaxFac* ...
              obj.simParameters.pTarget)];
          end
        end
        obj.dt = min([obj.dt * min(tmpVec),obj.simParameters.dtMax]);
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

    function getSystemSize(obj)
      % assemble blocks of jacobian matrix for multidomain system
      N = 0;
      for i = 1:obj.nDom
        nf = numel(obj.domains(i).DoFManager.getFieldList());
        N = N + nf;
      end
      Nfld = N;
      for i = 1:obj.nInterf
        N = N + obj.interfaces{i}.nFld;
      end
      Ninterf = N - Nfld;
      obj.systSize = [N,Nfld,Ninterf];
    end

    function getNumField(obj)
      for i = 1:obj.nDom
        obj.nfldDom(i) = numel(obj.domains(i).Discretizer.fields);
      end
      for i =1:obj.nInterf
        obj.nfldInt(i) = obj.interfaces{i}.nFld;
      end

    end

    function J = assembleJacobian(obj)
      % assemble blocks of jacobian matrix for multidomain system
      %[N,Nf,Ni] = deal(obj.systSize(1),obj.systSize(2),obj.systSize(3));
      J = cell(obj.systSize(1));
      f = 0;
      % populate jacobian with inner domain blocks
      for iD = 1:obj.nDom
        discr = obj.domains(iD).Discretizer;
        J(f+1:f+obj.nfldDom(iD),f+1:f+obj.nfldDom(iD)) = ...
          discr.assembleJacobian();
        for iF = 1:obj.nfldDom
          for iI = discr.interfaceList
            pos = find(strcmp(obj.interfaces{iI}.physics,discr.fields(iF)));
            jj = obj.systSize(2)+obj.nfldInt(iI)-obj.nfldInt(1)+pos;
            [J{iF+f,jj},J{jj,iF+f}] = getJacobian(...
              obj.interfaces{iF},pos,iD);
          end
        end
        f = f + obj.nfldDom(iD);  % update field counter
      end
    end

    function rhs = assembleRhs(obj)
      % assemble blocks of rhs for multidomain system
      rhs = cell(obj.systSize(1),1);
      f = 0;

      for iD = 1:obj.nDom
        discr = obj.domains(iD).Discretizer;
        rhs{f+1:f+obj.nfldDom(iD)} = discr.assembleRhs();
        for iF = 1:obj.nfldDom(iD)
          for iI = discr.interfaceList
            pos = find(strcmp(obj.interfaces{iI}.physics,discr.fields(iF)));
            rhs{f+iF} = rhs{f+iF}{:} + getRhs(...
              obj.interfaces{iF},pos,iD);
            iMult = obj.systSize(2)+sum(obj.nfldInt(2:iI))+pos;
            if isempty(rhs{iMult})
              % dont compute rhsMult twice: 1field -> 1 interface
              rhs{iMult} = getRhs(obj.interfaces{iF},pos);
            end        
          end
          f = f + 1;
        end
      end

    end
  end
end