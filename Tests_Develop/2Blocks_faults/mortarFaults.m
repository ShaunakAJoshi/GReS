classdef mortarFaults < handle
   
   properties
      simParameters
      t = 0
      nDom = 2
      nInt = 1
      tStep = 0
      iter
      dt
      state
      models
      meshGlue
   end
   
   methods
      function obj = mortarFaults(simParam,meshGlue)
         obj.setParams(simParam,meshGlue);
         obj.simulationLoop();
      end

      function simulationLoop(obj)
         obj.dt = obj.simParameters.dtIni;
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
         assert(numel(obj.models)==2,'Too many input domains');
         % state structure easier to access
         obj.state = repmat(struct('prev',{},'curr',{}),2,1);
         obj.state(1).prev = obj.meshGlue.model(1).State;
         obj.state(1).curr = copy(obj.state(1).prev);
         obj.state(2).prev = obj.meshGlue.model(1).State;
         obj.state(2).curr = copy(obj.state(1).prev);
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



   end
end

