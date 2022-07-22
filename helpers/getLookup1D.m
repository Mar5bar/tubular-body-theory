function lookup = getLookup1D()
%% Return a struct containing the indices in which the relevant fields can be
% found in the mesh array for a 1D mesh.
    lookup = struct();
    lookup.XBoundLower = 1;
    lookup.XBoundUpper = 2;
    lookup.XMid = 3;
end