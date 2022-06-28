function [mask, inds] = refinementNeeded(vals, neighbours, thresh)
    if nargin < 3
        thresh = 0.25;
    end
    numCells = length(neighbours);
    mask = zeros(numCells,1);
    for cellInd = 1 : numCells
        meanNeighbourVals = mean(vals(neighbours{cellInd}));
        mask(cellInd) = abs(vals(cellInd) - meanNeighbourVals) > thresh * meanNeighbourVals ;
    end

    if nargout == 2
        inds = find(mask);
    end
end