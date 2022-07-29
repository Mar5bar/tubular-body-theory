function drawCells(mesh2D)
    cells = mesh2D.cells;
    lookup = mesh2D.lookup;
    X = [cells(lookup.XBoundLower,:);cells(lookup.XBoundLower,:);cells(lookup.XBoundUpper,:);cells(lookup.XBoundUpper,:)];
    Y = [cells(lookup.YBoundLower,:);cells(lookup.YBoundUpper,:);cells(lookup.YBoundUpper,:);cells(lookup.YBoundLower,:)];
    patch(X,Y,'black','FaceColor','none')
end