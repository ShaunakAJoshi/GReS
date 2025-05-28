classdef SinglePhysics < handle
   properties (Access = public)
      J
      rhs
   end
   properties (SetAccess = private)
      field
   end
   properties
      model
      simParams
      dofm
      mesh
      elements
      faces
      material
      GaussPts
      state
   end
   
   methods
      function obj = SinglePhysics(fld,symmod,params,dofManager,grid,mat,state,data)
         obj.field = fld;
         obj.model = symmod;
         obj.simParams = params;
         obj.dofm = dofManager;
         obj.mesh = grid.topology;
         obj.elements = grid.cells;
         obj.faces = grid.faces;
         obj.material = mat;
         obj.state = state;
         % obj.field = 
         if ~isempty(data)
            obj.GaussPts = data{1};
         end
      end

      function applyNeuBC(obj,dofs,vals)
         % Base Neumann BCs application method
         % Vals are subtracted since we solve du = J\(-rhs)!
         obj.rhs(dofs) = obj.rhs(dofs) - vals; 
      end

      function applyDirBC(obj,~,dofs,varargin)
         % Base application of boundary condition to jacobian block. This
         % method only implements standard Dirichlet BCs on diagonal
         % jacobian blocks, but can be overridden by specific
         % implementations.
         % This version works with incremental linear system (vals = 0).
         % BC values for Dirichlet BC are zero (since the
         % system is solved in incremental form)
         % set Dir rows to zero
         obj.J = obj.J';
         obj.J(:,dofs) = 0; % setting columns is much faster
         obj.J = obj.J';
         % Update rhs with columns to be removed
         %obj.rhs = obj.rhs - obj.J(:,dofs)*vals;
         % set obj.rhs vals to dir vals
         obj.rhs(dofs) = 0;
         obj.J(:,dofs) = 0;
         % modify diagonal entries of K (avoiding for loop and accessing)
         Jdiag = diag(obj.J);
         Jdiag(dofs) = 1;
         obj.J = obj.J - diag(diag(obj.J)) + diag(Jdiag);
      end

      function J = getJacobian(obj,varargin)
         J = obj.J;
      end

      function rhs = getRhs(obj,varargin)
         rhs = obj.rhs;
      end


      function varargout = assembleFEM(obj,nEntries,assembler,N)
        % general sparse assembly loop over elements for single physics
        % kernel using FEM this method is flexible and can generate
        % multiple sparse matrices in output. nEntries: cell array of
        % length of index arrays for sparse assembly 
        % localAssembly: cell array of routines to compute local 
        % N: cell array of [numb row,numb col] for each output marix

        % setup indices array
        assert(numel(nEntries)==numel(assembler),...
          'Lenght of entires number and assembly routine list must match!');

        subCells = obj.dofm.getFieldCells(obj.field);
        nSubCellsByType = histc(obj.mesh.cellVTKType(subCells),[10, 12, 13, 14]);
        % number of matrices to be assembled
        nMat = numel(nEntries);
        assert(nargout == nMat,'Incorrect number of output matrices.')
        varargout = cell(nMat,1);
        [iiVec,jjVec,matVec] = deal(cell(nMat,1));

        for i = 1:nMat
          % preallocate arrays for sparse assembly
          n = nEntries{i}*nSubCellsByType;
          [iiVec{i},jjVec{i},matVec{i}] = deal(zeros(n,1));
        end
        l1 = 0;
        l2 = 0;
        % loop over cells
        for el = subCells'
          for i = 1:nMat
            % get dof id and local matrix
            [dofRow,dofCol,locMat,s2] = assembler{i}(el,l2);
            [jjLoc,iiLoc] = meshgrid(dofCol,dofRow);
            s1 = numel(locMat(:));
            iiVec{i}(l1+1:l1+s1) = iiLoc(:);
            jjVec{i}(l1+1:l1+s1) = jjLoc(:);
            matVec{i}(l1+1:l1+s1) = locMat(:);
            l1 = l1 + s1;
            l2 = l2 + s2;
          end
        end

        for i = 1:nMat
          % renumber indices according to active nodes
          % important: this call to unique assumes that iiVec contains all active
          % degrees of freedom in the domain
          [~,~,iiVec{i}] = unique(iiVec{i});
          [~,~,jjVec{i}] = unique(jjVec{i});
          % populate stiffness matrix
          Nr = N{i}(1); 
          Nc = N{i}(2);
          varargout{i} = sparse(iiVec{i}, jjVec{i}, matVec{i}, Nr, Nc);
        end
      end
   end
end