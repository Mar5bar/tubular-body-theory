function M = symmMat(entries,dim)
%% symmMat constructs a symmetric matrix of size [dim,dim] with lower
%% triangular entries given in linear order, eg entries = [M11, M12, M22] for a 2x2 matrix.
    assert(numel(entries) == dim*(dim+1)/2, 'Number of entries incompatible with matrix size.');
    M = zeros(dim);
    diagInds = cumsum([1,dim:-1:2]);
    for col = 1 : dim - 1
        M(col+1:dim,col) = entries(diagInds(col) + 1 : diagInds(col + 1) - 1);
    end
    M = M + transpose(M) + diag(entries(diagInds));
end