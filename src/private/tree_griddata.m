function [xx,yy,vv] = tree_griddata(tree)
xx = []; yy = []; vv = [];
cleaves = tree.leaves();
for lvcnt = 1:length(cleaves)
    cleaf = cleaves{lvcnt};
    % GRID POINTS
    resPerNode = cleaf.data.resolution;
    [xr,yr,zr,dx,dy,dz] = cleaf.mesh(resPerNode);
    vals = cleaf.data.values;
    xx = [xx; xr(:)];
    yy = [yy; yr(:)];
    vv = [vv; vals(:)];
end
end