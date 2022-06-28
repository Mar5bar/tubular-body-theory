function [cells, numCells, neighbours] = refineCells(cells, lookup, numCells, refInds)

    % For each cell in refInds, split the cell into four pieces and delete the
    % original cell.
    newCell = zeros(6,1);
    for i = 1 : numel(refInds)
        ind = refInds(i);
        oldCell = cells(:,ind);

        % Create the upper-left cell.
        numCells = numCells + 1;

        newCell(lookup.XBoundLower) = oldCell(lookup.XBoundLower);
        newCell(lookup.XBoundUpper) = oldCell(lookup.XMid);
        newCell(lookup.YBoundLower) = oldCell(lookup.YMid);
        newCell(lookup.YBoundUpper) = oldCell(lookup.YBoundUpper);
        newCell(lookup.XMid) = mean([newCell(lookup.XBoundLower), newCell(lookup.XBoundUpper)]);
        newCell(lookup.YMid) = mean([newCell(lookup.YBoundLower), newCell(lookup.YBoundUpper)]);

        cells(:,numCells) = newCell;

        % Create the upper-right cell.
        numCells = numCells + 1;

        newCell(lookup.XBoundLower) = oldCell(lookup.XMid);
        newCell(lookup.XBoundUpper) = oldCell(lookup.XBoundUpper);
        newCell(lookup.YBoundLower) = oldCell(lookup.YMid);
        newCell(lookup.YBoundUpper) = oldCell(lookup.YBoundUpper);
        newCell(lookup.XMid) = mean([newCell(lookup.XBoundLower), newCell(lookup.XBoundUpper)]);
        newCell(lookup.YMid) = mean([newCell(lookup.YBoundLower), newCell(lookup.YBoundUpper)]);

        cells(:,numCells) = newCell;

        % Create the lower-left cell.
        numCells = numCells + 1;

        newCell(lookup.XBoundLower) = oldCell(lookup.XBoundLower);
        newCell(lookup.XBoundUpper) = oldCell(lookup.XMid);
        newCell(lookup.YBoundLower) = oldCell(lookup.YBoundLower);
        newCell(lookup.YBoundUpper) = oldCell(lookup.YMid);
        newCell(lookup.XMid) = mean([newCell(lookup.XBoundLower), newCell(lookup.XBoundUpper)]);
        newCell(lookup.YMid) = mean([newCell(lookup.YBoundLower), newCell(lookup.YBoundUpper)]);

        cells(:,numCells) = newCell;

        % Create the lower-right cell.
        numCells = numCells + 1;

        newCell(lookup.XBoundLower) = oldCell(lookup.XMid);
        newCell(lookup.XBoundUpper) = oldCell(lookup.XBoundUpper);
        newCell(lookup.YBoundLower) = oldCell(lookup.YBoundLower);
        newCell(lookup.YBoundUpper) = oldCell(lookup.YMid);
        newCell(lookup.XMid) = mean([newCell(lookup.XBoundLower), newCell(lookup.XBoundUpper)]);
        newCell(lookup.YMid) = mean([newCell(lookup.YBoundLower), newCell(lookup.YBoundUpper)]);

        cells(:,numCells) = newCell;

    end

    % Delete the old cells.
    cells(:,refInds) = [];
    numCells = numCells - numel(refInds);

    % Compute the new neighbours.
    neighbours = genNeighbours(cells, lookup, numCells);
end