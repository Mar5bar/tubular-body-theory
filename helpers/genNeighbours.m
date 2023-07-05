function neighbours = genNeighbours(cells, lookup, numCells)
% Generate a list of all neighbours.
    neighbours = cell(numCells,1);
    for cellInd = 1 : numCells
        myXBoundLower = cells(lookup.XBoundLower,cellInd);
        myXBoundUpper = cells(lookup.XBoundUpper,cellInd);
        myYBoundLower = cells(lookup.YBoundLower,cellInd);
        myYBoundUpper = cells(lookup.YBoundUpper,cellInd);
        
        xMask = isWithin(cells(lookup.XBoundLower,1:numCells), myXBoundLower, myXBoundUpper) | isWithin(cells(lookup.XBoundUpper,1:numCells), myXBoundLower, myXBoundUpper);
        yMask = isWithin(cells(lookup.YBoundLower,1:numCells), myYBoundLower, myYBoundUpper) | isWithin(cells(lookup.YBoundUpper,1:numCells), myYBoundLower, myYBoundUpper);
        neighbours{cellInd} = find(xMask & yMask);
    end
end
