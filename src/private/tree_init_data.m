function tree_init_data(tree, func, resPerNode)
cleaves = tree.leaves();
for lvcnt = 1:length(cleaves)
    cleaf = cleaves{lvcnt};
    [xr,yr,zr,dx,dy,dz] = cleaf.mesh(resPerNode);
    cvalues = func(xr,yr);
    % TODO: get the correct dimension
    cleaf.data.dim    = 1;
    cleaf.data.resolution = resPerNode;
    cleaf.data.values = cvalues;
end
end