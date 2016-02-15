%/* ************************************************** */
%     Copyright 2014 Arash Bakhtiari
%
%     you may not use this file except in compliance with the License.
%     You obtain a copy of the License in the LICENSE file
%
%     Unless required by applicable law or agreed to in writing, software
%     distributed under the License is distributed on an "AS IS"" BASIS,
%     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%     See the License for the specific language governing permissions and
%     limitations under the License.
%/* ************************************************** */

classdef tree_data < handle
%UNTITLED Summary of this class goes here
%   Detailed explanation goes here

    properties
    end

    methods (Static)

        %/* ************************************************** */
        function init_data(tree, func, resPerNode, t)
        % Initializes the values of all leaf octants based on the
        % input function values and grid resolution per node
            if nargin < 4, t = 0; end;
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                tree_data.set_node_fn(cleaf, func, resPerNode, t);
            end
        end

        %/* ************************************************** */
        function interp_tree(src_tree, dst_tree)
        % interpolate the values of grid points in the destination tree
        % from the values of correspoding nodes in the source tree
        % TODO: this is a brute-force algorithm. optimize it.
            global verbose;
            global INTERP_TYPE;
            global resPerNode;
            dst_leaves  = dst_tree.leaves();
            for dst_lvcnt = 1:length(dst_leaves)
                dst_leaf = dst_leaves{dst_lvcnt};
                [xx,yy,zz,dx,dy,dz] = dst_leaf.mesh(resPerNode, INTERP_TYPE);
                tmpval = tree_data.interp_points(src_tree, xx, yy, zz, INTERP_TYPE);
                tree_data.set_node_val(dst_leaf, tmpval, resPerNode);
            end
        end

        %/* ************************************************** */
        function val = interp_points(src_tree,xq,yq,zq,INTERP_TYPE)
            src_leaves  = src_tree.leaves();
            resPerNode  = src_leaves{1}.data.resolution;
            data_dim    = src_leaves{1}.data.dim;
            val = zeros([size(xq) 1 data_dim]);
            for dimcnt = 1:data_dim
                % interpolate the destination values from corresponding source leaves
                tmpval = zeros(size(xq));
                for src_lvcnt =1:length(src_leaves)
                    src_leaf = src_leaves{src_lvcnt};

                    % find the points inside current source leaf
                    indices = tree_data.points_in_node(src_leaf, xq, yq);
                    if ~any(any(indices)), continue; end;

                    xx = xq(indices);
                    yy = yq(indices);

                    if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                        global CHEB_IMPL
                        if strcmp(CHEB_IMPL, 'IAS')
                            vv = zeros(size(xx));
                            w = src_leaf.data.values;
                            [xmin,xmax,ymin,ymax] = ...
                                src_leaf.corners;
                            xs = (xx - xmin)*2/(xmax-xmin)-1.0;
                            ys = (yy - ymin)*2/(ymax-ymin)-1.0;
                            for xindx =1:size(xx)
                                vv(xindx) =  tree_data.chebeval2(w, xs(xindx), ys(xindx));
                            end
                        elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                            w = src_leaf.data.values;
                            vv = w(xx,yy);
                        end
                    else
                        [src_xxr,src_yyr,src_zzr,src_dx,src_dy,src_dz] = src_leaf.mesh(resPerNode);
                        interp_data = src_leaf.data.values(:,:,:,dimcnt);
                        vv = interp2(src_xxr,src_yyr,interp_data,xx,yy, ...
                                     INTERP_TYPE);
                    end
                    tmpval(indices) = vv;
                end
                val(:,:,:,dimcnt) = tmpval;
            end
        end

        %/* ************************************************** */
        function max_err = compute_error(tree, fexact, t,INTERP_TYPE)
            leaves     = tree.leaves();
            resPerNode = leaves{1}.data.resolution;
            data_dim   = leaves{1}.data.dim;
            max_err = 0;
            max_err_pnt = [];
            max_err_node = [];
            for lvcnt =1:length(leaves)
                leaf = leaves{lvcnt};
                [xr,yr,zr,dx,dy,dz] = leaf.mesh(resPerNode);
                % compute the center of the local grid cells
                xxc = xr+dx/2;
                yyc = yr+dy/2;
                zzc = zr+dz/2;
                xxcc = xxc(1:end-1,1:end-1);
                yycc = yyc(1:end-1,1:end-1);
                zzcc = zzc(1:end-1,1:end-1);

                vale = fexact(t,xxcc,yycc,zzcc);
                valt = tree_data.interp_points(leaf,xxcc,yycc,zzcc,INTERP_TYPE);

                diff = vale - valt;
                [err, indx] = max(abs(diff(:)));
                [i_row, i_col] = ind2sub(size(diff),indx);
                if err > max_err
                    max_err = err;
                    max_err_node = leaf;
                    max_err_pnt = [i_row, i_col];
                end;
            end
        end

        % %/* ************************************************** */
        % function [tree_out] = collapse(tree_in)
        % % TODO: check that treecells have the exact same structure
        % %       -> mids of all treecells are the same.
        % % clone structure of the resulting tree from the input trees
        %     num_trees = length(tree_in);
        %     num_leaves = length(tree_in{1}.leaves());
        %     tree_out = qtree.clone(tree_in{1});

        %     % get the leaves of input trees
        %     tree_in_leaves = cell(num_trees,num_leaves);
        %     for tree_in_cnt =1:num_trees
        %         tree_in_leaves(tree_in_cnt,:) = tree_in{tree_in_cnt}.leaves();
        %     end

        %     % iterate over leaves of tree out
        %     tree_out_leaves = tree_out.leaves();
        %     for leafcnt = 1:length(tree_out_leaves)
        %         leaf = tree_out_leaves{leafcnt};

        %         tmp = tree_in_leaves{1,leafcnt};
        %         leaf.data.dim           = num_trees;
        %         leaf.data.resolution    = tmp.data.resolution;
        %         % TODO: remove the one after extending the code to 3D
        %         leaf.data.values        = zeros([size(tmp.data.values) 1 num_trees]);

        %         for tree_in_cnt = 1:num_trees
        %             tree_in_leaf = tree_in_leaves{tree_in_cnt,leafcnt};
        %             leaf.data.values(:,:,:,tree_in_cnt) = tree_in_leaf.data.values;
        %         end
        %     end
        % end

        %/* ************************************************** */
        function set_node_fn(node, fn, resPerNode, t)
            global INTERP_TYPE
            [xx,yy,zz,dx,dy,dz] = node.mesh(resPerNode, INTERP_TYPE);
            fnval = fn(t,xx,yy,zz);
            tree_data.set_node_val(node, fnval, resPerNode);
        end

        %/* ************************************************** */
        function set_node_val(node, fnval, resPerNode)
            global INTERP_TYPE
            if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                w = tree_data.get_node_cheb_interpolant(node, fnval, resPerNode);
                node.data.values = w;
            else
                node.data.values = fnval;
            end
            node.data.dim = 1;
            node.data.resolution = resPerNode;
        end

        %/* ************************************************** */
        function [w] = get_node_cheb_interpolant(node, fn_val, ...
                                                       resPerNode)
            n1 = resPerNode+1;
            n2 = resPerNode+1;
            global CHEB_IMPL;
            if strcmp(CHEB_IMPL, 'IAS')
                % roots grid
                i = [0:(n1-1)]';
                j = [0:(n2-1)]';
                x = cos((i+1/2)*pi/n1);
                y = cos((j+1/2)*pi/n2);
                % compute the value of chebyshev polynomials at the node points.
                Tx = tree_data.chebpoly(n1-1,x);
                Ty = tree_data.chebpoly(n2-1,y);
                % using discrete orthogonality of chebyshev polynomials to compute the
                % coefficients for approximation using chebyshev polynomial basis.
                for i=1:n1
                    w_(i,:) = (reshape(fn_val(i,:),1,[])*Ty)*2/n2;
                end
                for j=1:n2
                    w(:,j) = (reshape(w_(:,j),1,[])*Tx)*2/n1;
                end
                w(1,:) = w(1,:)/2;
                w(:,1) = w(:,1)/2;
            elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                [xmin xmax ymin ymax] = node.corners;
                w = chebfun2(fn_val, [xmin xmax ymin ymax]);
            end
        end

        %/* ************************************************** */
        function [val] = get_node_values_cheb(node, resPerNode)
            n1 = resPerNode+1;
            n2 = resPerNode+1;
            global CHEB_IMPL;
            if strcmp(CHEB_IMPL, 'IAS')
                % roots grid
                i = [0:(n1-1)]';
                j = [0:(n2-1)]';
                x = cos((i+1/2)*pi/n1);
                y = cos((j+1/2)*pi/n2);
                w = node.data.values;
                val = tree_data.chebeval2(w,x,y);
                % CHEBFUN IMPLEMENTATION
            elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                w = node.data.values;
                [xx,yy,zz,dx,dy,dz] = node.mesh(resPerNode, ...
                                                'CHEBYSHEV');
                val = w(xx,yy);
            end
        end

        %/* ************************************************** */
        function [f_val] = chebeval2(w,x,y)
        %Compute the values of a chebyshev approximation at a
        %regular grid specified by x,y.
            f_val = zeros(length(x),length(y));
            T_x = tree_data.chebpoly(size(w,1)-1, x);
            T_y = tree_data.chebpoly(size(w,2)-1, y);
            for i=1:size(w,1)
                f_(i,:)=T_y*reshape(w(i,:),[],1);
            end
            for j=1:length(y)
                f_val(:,j)=T_x*reshape(f_(:,j),[],1);
            end
        end

        %/* ************************************************** */
        function [T] = chebpoly(n,x)
        %Retruns values of all chebyshev polynomials upto degree n
        %at points x
            x_ = reshape(x,[],1);
            T0 = ones(size(x_));
            T(:,1) = T0;
            if n == 0
                return;
            end
            T1 = x_;
            T(:,2) = T1;
            if n == 1
                return;
            end
            for i = 2:n
                T(:,i+1) = 2*x_.*T1-T0;
                T0 = T1;
                T1 = T(:,i+1);
            end
        end

        %/* ************************************************** */
        function [xx,yy,vv] = grid_points(tree)
            global INTERP_TYPE;
            global resPerNode;
            xx = []; yy = [];
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                if isempty(cleaf.data), continue; end;
                [xr,yr,zr,dx,dy,dz] = cleaf.mesh(resPerNode,INTERP_TYPE);
                xx = [xx; xr(:)];
                yy = [yy; yr(:)];
            end
        end

        %/* ************************************************** */
        function [vv] = grid_data(tree)
            global INTERP_TYPE;
            vv = [];
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                if isempty(cleaf.data), continue; end;
                % GRID VALUES
                if strcmp(INTERP_TYPE,'CHEBYSHEV')
                    %evaluate chebyshev polynomials at the cheb
                    resPerNode = cleaf.data.resolution;
                    vals = tree_data.get_node_values_cheb(cleaf, resPerNode);
                    tmp = vals(:,:);
                else
                    vals = cleaf.data.values;
                    tmp = vals(:,:,:);
                end
                vv = [vv; tmp(:)];
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
            %tree.plottree(0.5)
            hold on
            [txx,tyy]   = tree_data.grid_points(tree);
            [tvv]       = tree_data.grid_data(tree);
            if ~isempty(tvv) & ~isempty(txx) & ~isempty(tyy)
                scatter3(txx,tyy,tvv(:,dim),ones(size(txx)),tvv(:,dim),'filled')
            end
            %plot3(txx(:),tyy(:),tvv(:,dim),'.',MS,1)%,ones(size(txx)),tvv(:,dim),'filled')
            axis off;
            axis equal;
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
end
