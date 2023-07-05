function lookup = getLookup2D()
%% Return a struct containing the indices in which the relevant fields can be
% found in the mesh array for a 2D mesh.
    lookup = struct();
    lookup.XBoundLower = 1;
    lookup.XBoundUpper = 2;
    lookup.YBoundLower = 3;
    lookup.YBoundUpper = 4;
    lookup.XMid = 5;
    lookup.YMid = 6;
    lookup.XWidth = 7;
    lookup.YWidth = 8;
end