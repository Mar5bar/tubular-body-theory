function [cells, lookup, numCells, neighbours] = initMesh(bounds, initXNum, initYNum, maxCells)
    cells = zeros(6,maxCells);
    % Lookup for indices in each cell.
    lookup = struct();
    lookup.XBoundLower = 1;
    lookup.XBoundUpper = 2;
    lookup.YBoundLower = 3;
    lookup.YBoundUpper = 4;
    lookup.XMid = 5;
    lookup.YMid = 6;

    % Form the initial discretisation.
    XEnds = linspace(bounds.X(1), bounds.X(2), initXNum+1);
    YEnds = linspace(bounds.Y(1), bounds.Y(2), initYNum+1);
    XMids = movmean(XEnds,2,'Endpoints','discard');
    YMids = movmean(YEnds,2,'Endpoints','discard');

    numCells = 0;
    for Xind = 1 : initXNum
        for Yind = 1 : initYNum
            numCells = numCells + 1;
            
            newCell = zeros(6,1);

            newCell(lookup.XBoundLower) = XEnds(Xind);
            newCell(lookup.XBoundUpper) = XEnds(Xind+1);
            newCell(lookup.YBoundLower) = YEnds(Yind);
            newCell(lookup.YBoundUpper) = YEnds(Yind+1);
            newCell(lookup.XMid) = XMids(Xind);
            newCell(lookup.YMid) = YMids(Yind);

            cells(:,numCells) = newCell;
        end
    end
    % Generate a list of neighbours for each cell.
    neighbours = genNeighbours(cells, lookup, numCells);
end