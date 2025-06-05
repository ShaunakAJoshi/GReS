classdef Elements < handle
  % ELEMENTS General element class

  properties (Access = public)
    % Mesh reference
    mesh

    % Element and surface type counts
    nCellsByType     % [#tetra, #hexa, #wed, #pyr]
    nSurfByType      % [#tri, #quad]
    nNodesElem = [4, 8, 6, 5]

    % Geometric properties
    centroid
    coordLoc
    
    vol

    % 3D element handlers
    tetra
    hexa

    % 2D element handlers
    tri
    quad
    quadL

    % Derivative indices
    indB
    indB2D
  end

  properties (Access = private)
    GaussPts = []
  end

  methods (Access = public)
    function obj = Elements(mesh, varargin)
      data = varargin;
      obj.setElementData(nargin, mesh, data);
      obj.computeCellProperties;
    end

    function computeCellProperties(obj)
      idCell = 1:obj.mesh.nCells;
      idTetra = idCell(obj.mesh.cellVTKType == 10);
      idHexa  = idCell(obj.mesh.cellVTKType == 12);

      if ~isempty(idTetra)
        obj.vol(idTetra) = obj.tetra.findVolume(idTetra);
        obj.cellCentroid(idTetra,:) = computeCentroidGeneral(obj, idTetra);
      end

      if ~isempty(idHexa)
        [obj.vol(idHexa), obj.cellCentroid(idHexa,:)] = ...
          obj.hexa.findVolumeAndCentroid(idHexa);
      end
    end

    function elem = createElement(obj, id)
      switch id
        case 'Hexahedron'
          elem = Hexahedron(obj.mesh, obj.GaussPts);
        case 'Tetrahedron'
          elem = Tetrahedron(obj.mesh, obj.GaussPts);
        case 'Quadrilateral'
          elem = Quadrilateral(obj.mesh, obj.GaussPts);
        case 'Triangle'
          elem = Triangle(obj.mesh, obj.GaussPts);
      end
    end
  end

  methods (Access = private)
    function setElementData(obj, nIn, mesh, data)
      obj.mesh = mesh;
      if nIn > 1
        obj.GaussPts = data{1};
      end

      % Count 3D element types
      obj.nCellsByType = histc(obj.mesh.cellVTKType, [10, 12, 13, 14]);
      if isempty(obj.nCellsByType)
        obj.nCellsByType = zeros(4,1);
      end

      % Instantiate 3D element classes
      if obj.nCellsByType(1) > 0
        obj.tetra = Tetrahedron(obj.mesh);
      end
      if obj.nCellsByType(2) > 0
        obj.hexa = Hexahedron(obj.mesh, obj.GaussPts);
      end

      obj.vol = zeros(obj.mesh.nCells, 1);
      obj.cellCentroid = zeros(obj.mesh.nCells, 3);

      % Count 2D surface types
      obj.nSurfByType = histc(obj.mesh.surfaceVTKType, [5, 9]);

      if obj.nSurfByType(1) > 0
        obj.tri = Triangle(obj.mesh, obj.GaussPts);
      end

      if obj.nSurfByType(2) > 0
        if any(obj.mesh.surfaceNumVerts == 4)
          obj.quad = Quadrilateral(obj.mesh, obj.GaussPts);
        elseif any(obj.mesh.surfaceNumVerts == 8)
          obj.quad = Quad8(obj.mesh, obj.GaussPts);
          if obj.mesh.cartGrid
            obj.quadL = Quadrilateral(getQuad4mesh(obj.mesh), obj.GaussPts);
          end
        end
      end

      % Build indB2D
      if obj.nSurfByType(2) == 0
        l1 = 3;
      else
        l1 = 4 * obj.GaussPts.nNode;
      end

      obj.indB2D = zeros(4 * l1, 2);
      obj.indB2D(:,1) = repmat([1, 2, 2, 1], [1, l1]);
      obj.indB2D(:,2) = repmat([1, 3, 5, 6], [1, l1]);
      obj.indB2D(:,1) = obj.indB2D(:,1) + repelem(2*(0:(l1-1))', 4);
      obj.indB2D(:,2) = obj.indB2D(:,2) + repelem(6*(0:(l1-1))', 4);
    end

    function tetraCtr = computeCentroidGeneral(obj, idTetra)
      tetraCtr = sparse(...
        repelem(1:length(idTetra), obj.mesh.cellNumVerts(idTetra)), ...
        nonzeros(obj.mesh.cells(idTetra,:)'), ...
        repelem(1 ./ obj.mesh.cellNumVerts(idTetra), obj.mesh.cellNumVerts(idTetra)), ...
        length(idTetra), obj.mesh.nNodes) * obj.mesh.coordinates;
    end
  end

  methods (Static)

    function centroidCoord = getElemCentroid(vtkType)
      switch vtkType
        case 5
          centroidCoord = Triangle.centroid;
        case 9
          centroidCoord = Quadrilateral.centroid;
        case 12
          centroidCoord = Hexahedron.centroid;
        otherwise
          error('centroid retrieving for element vtk type %i not yet implemented',vtkType)
      end
    end

  end
end
