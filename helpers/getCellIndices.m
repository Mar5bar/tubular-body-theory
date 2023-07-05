function [cellInd1, cellInd2] = getCellIndices(l, numCells)
%% getCellIndices takes a linear index of numCells^2 ordered pairs of cells
%% and returns the indices of the two cells.
    [cellInd1, cellInd2] = ind2sub([numCells, numCells], l);
end