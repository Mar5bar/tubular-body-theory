function drawCells(mesh2D,color)
    if nargin < 2
        color = 'black';
    end
    cells = mesh2D.cells(:,1:mesh2D.numCells);
    lookup = mesh2D.lookup;
    X = [cells(lookup.XBoundLower,:);cells(lookup.XBoundLower,:);cells(lookup.XBoundUpper,:);cells(lookup.XBoundUpper,:)];
    Y = [cells(lookup.YBoundLower,:);cells(lookup.YBoundUpper,:);cells(lookup.YBoundUpper,:);cells(lookup.YBoundLower,:)];
    patch(X,Y,'black','FaceColor','none','EdgeColor',color)
end