classdef tree_data < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        %/* ************************************************** */
        function [xx,yy,vv] = grid_points(tree)
            xx = []; yy = [];
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                % GRID POINTS
                resPerNode = cleaf.data.resolution;
                [xr,yr,zr,dx,dy,dz] = cleaf.mesh(resPerNode);
                xx = [xx; xr(:)];
                yy = [yy; yr(:)];
            end
        end
        
        %/* ************************************************** */
        function [vv] = grid_data(tree)
            vv = [];
            cleaves = tree.leaves();
            dim = cleaves{1}.data.dim;
            for dcnt = 1:dim
                vvtmp = [];
                for lvcnt = 1:length(cleaves)
                    cleaf = cleaves{lvcnt};                 
                    % GRID VALUES                    
                    vals = cleaf.data.values;                                       
                    tmp = vals(:,:,:,dcnt);
                    vvtmp = [vvtmp; tmp(:)];
                end
                vv(:,dcnt) = vvtmp;
            end
        end
        
        %/* ************************************************** */
        function init_data(tree, func, resPerNode)
            
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                [xr,yr,zr,dx,dy,dz] = cleaf.mesh(resPerNode);
                cvalues = func(xr,yr);
                cleaf.data.dim          = size(cvalues,4);
                cleaf.data.resolution   = resPerNode;
                for dcnt = 1:size(cvalues,4)
                    cleaf.data.values(:,:,:,dcnt)   = cvalues(:,:,:,dcnt);
                end
            end
        end
        
        %/* ************************************************** */
        function interp(src_tree, dst_tree)
            % interpolate the values of grid points in the destination tree
            % from the values of correspoding nodes in the source tree
            % TODO: this is a brute-force algorithm. optimize it.            
            global verbose;
            
            dst_leaves  = dst_tree.leaves();
            src_leaves  = src_tree.leaves();
            
            resPerNode  = src_leaves{1}.data.resolution;
            data_dim    = src_leaves{1}.data.dim;
            
            for dst_lvcnt = 1:length(dst_leaves)
                dst_leaf = dst_leaves{dst_lvcnt};
                
                % get the mesh of destination leaf
                [dst_xxr,dst_yyr,dst_zzr,dst_dx,dst_dy,dst_dz] = dst_leaf.mesh(resPerNode);
                
                % init the destination leaf data values
                dst_leaf.data.dim = data_dim;
                dst_leaf.data.resolution = resPerNode;
                % remove the 1 for 3D implementation
                dst_leaf.data.values = zeros([size(dst_xxr) 1 data_dim]);
                
                for dimcnt = 1:data_dim
                    % interpolate the destination values from corresponding source leaves
                    tmpval = zeros(size(dst_xxr));
                    for src_lvcnt =1:length(src_leaves)
                        src_leaf = src_leaves{src_lvcnt};
                        
                        % find the points inside current source leaf
                        indices = points_in_node(src_leaf, dst_xxr, dst_yyr);
                        if ~any(any(indices)), continue; end;
                        
                        xx = dst_xxr(indices);
                        yy = dst_yyr(indices);
                        if verbose
                            fprintf('interpolate values of dst node %d from src node %d\n',dst_lvcnt, src_lvcnt);
                        end
                                                
                        [src_xxr,src_yyr,src_zzr,src_dx,src_dy,src_dz] = src_leaf.mesh(resPerNode);
                                                
                        interp_data = src_leaf.data.values(:,:,:,dimcnt);
                        vv = interp2(src_xxr,src_yyr,interp_data,xx,yy);
                        tmpval(indices) = vv;                                                
                    end
                    dst_leaf.data.values(:,:,:,dimcnt) = tmpval;
                end
            end
            
            %/* ************************************************** */
            function indices = points_in_node(node, xx, yy)
                % complication for points that lie right on
                % the boundaries
                [xmin,xmax,ymin,ymax]=corners(node);
                idx = find(xmin <= xx & xx <= xmax);
                idy = find(ymin <= yy & yy <= ymax);
                indices = intersect(idx, idy);
            end
        end
        
        %/* ************************************************** */
        function [tree_out] = collapse(tree_in)
            % TODO: check that treecells have the exact same structure
            %       -> mids of all treecells are the same.
            % clone structure of the resulting tree from the input trees
            num_trees = length(tree_in);
            num_leaves = length(tree_in{1}.leaves());
            tree_out = qtree.clone(tree_in{1});
            
            % get the leaves of input trees
            tree_in_leaves = cell(num_trees,num_leaves);
            for tree_in_cnt =1:num_trees
                tree_in_leaves(tree_in_cnt,:) = tree_in{tree_in_cnt}.leaves();
            end
            
            % iterate over leaves of tree out
            tree_out_leaves = tree_out.leaves();
            for leafcnt = 1:length(tree_out_leaves)
                leaf = tree_out_leaves{leafcnt};
                
                tmp = tree_in_leaves{1,leafcnt};
                leaf.data.dim           = num_trees;
                leaf.data.resolution    = tmp.data.resolution;
                % TODO: remove the one after extending the code to 3D
                leaf.data.values        = zeros([size(tmp.data.values) 1 num_trees]);
                
                for tree_in_cnt = 1:num_trees
                    tree_in_leaf = tree_in_leaves{tree_in_cnt,leafcnt};
                    leaf.data.values(:,:,:,tree_in_cnt) = tree_in_leaf.data.values;
                end
            end            
        end
        
        %/* ************************************************** */
        function plot_grid(tree)
            MS='MarkerSize';
            [txx,tyy] = tree_data.grid_points(tree);
            tree.plottree;
            hold on;
            plot(txx(:),tyy(:),'ro',MS,1);
            axis off; axis equal;
        end
        
        %/* ************************************************** */
        function plot_data(tree,dim)
            if nargin < 2, dim = 1; end;
            MS='MarkerSize';
            tree.plottree
            hold on
            [txx,tyy]   = tree_data.grid_points(tree);
            [tvv]       = tree_data.grid_data(tree);
            scatter3(txx,tyy,tvv(:,dim),ones(size(txx)),tvv(:,dim),'filled')                        
            axis off; 
            axis equal;
        end
    end
end

