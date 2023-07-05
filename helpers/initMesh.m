function m = initMesh(XBounds, YBounds, initXNum, initYNum, maxCells)

    m = struct();
    m.XBounds = XBounds;
    m.YBounds = YBounds;
    m.cells = zeros(8,maxCells);
    % Load the indices needed for looking up particular quantities in each cell.
    m.lookup = getLookup2D();

    % Form the initial discretisation.
    XEnds = linspace(XBounds(1), XBounds(2), initXNum+1);
    YEnds = linspace(YBounds(1), YBounds(2), initYNum+1);
    XMids = movmean(XEnds,2,'Endpoints','discard');
    YMids = movmean(YEnds,2,'Endpoints','discard');

    numCells = 0;
    for Xind = 1 : initXNum
        for Yind = 1 : initYNum
            numCells = numCells + 1;
            
            newCell = zeros(8,1);

            newCell(m.lookup.XBoundLower) = XEnds(Xind);
            newCell(m.lookup.XBoundUpper) = XEnds(Xind+1);
            newCell(m.lookup.YBoundLower) = YEnds(Yind);
            newCell(m.lookup.YBoundUpper) = YEnds(Yind+1);
            newCell(m.lookup.XMid) = XMids(Xind);
            newCell(m.lookup.YMid) = YMids(Yind);
            newCell(m.lookup.XWidth) = newCell(m.lookup.XBoundUpper) - newCell(m.lookup.XBoundLower);
            newCell(m.lookup.YWidth) = newCell(m.lookup.YBoundUpper) - newCell(m.lookup.YBoundLower);

            m.cells(:,numCells) = newCell;
        end
    end
    m.numCells = numCells;

    % Generate a list of neighbours for each cell.
    m.neighbours = genNeighbours(m.cells, m.lookup, numCells);
end