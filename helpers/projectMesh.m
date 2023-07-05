function mesh1D = projectMesh(mesh2D)
%% projectMesh takes in a 2D mesh structure and flattens it onto a 1D mesh for
% the X coordinate. mesh1D includes a map from mesh2D onto the 1D cell(s) that
% make up the 2D cell.
    
    mesh1D = struct();
    lookup = getLookup1D();

    % Use the X boundaries of the 2D mesh to form a 1D mesh.
    XBounds = unique(mesh2D.cells([mesh2D.lookup.XBoundLower,mesh2D.lookup.XBoundUpper],1:mesh2D.numCells));
    numCells = numel(XBounds) - 1;
    mesh1D.cells = zeros(4, numCells);
    mesh1D.cells(lookup.XBoundLower,:) = XBounds(1:end-1);
    mesh1D.cells(lookup.XBoundUpper,:) = XBounds(2:end);
    mesh1D.cells(lookup.XMid,:) = movmean(XBounds,2,'Endpoints','discard');
    mesh1D.cells(lookup.XWidth,:) = mesh1D.cells(lookup.XBoundUpper,:) - mesh1D.cells(lookup.XBoundLower,:);

    mesh1D.lookup = lookup;
    mesh1D.numCells = numCells;

    % We'll need to be able to take a 2D cell and identify which 1D cell(s) it
    % is made up of. We'll do this in a projMap, which will take the index of
    % a 2D cell and give the indices of the 1D cells that it corresponds to.
    mesh1D.projMap = cell(mesh2D.numCells,1);
    for cellInd = 1 : mesh2D.numCells
        first1DCell = find(mesh1D.cells(lookup.XBoundLower,:) == mesh2D.cells(mesh2D.lookup.XBoundLower,cellInd),1,'first');
        last1DCell = find(mesh1D.cells(lookup.XBoundUpper,:) == mesh2D.cells(mesh2D.lookup.XBoundUpper,cellInd),1,'first');
        mesh1D.projMap{cellInd} = first1DCell : last1DCell;
    end

end