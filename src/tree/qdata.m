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

classdef qdata < handle
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
                qdata.set_node_fn(cleaf, func, resPerNode, t);
            end
        end

        %/* ************************************************** */
        function interp_tree(src_tree, dst_tree)
        % interpolate the values of grid points in the destination tree
        % from the values of correspoding nodes in the source tree
            global VERBOSE;
            global INTERP_TYPE;
            global RES_PER_NODE;
            dst_leaves  = dst_tree.leaves();
            for dst_lvcnt = 1:length(dst_leaves)
                dst_leaf = dst_leaves{dst_lvcnt};
                [xx,yy,zz,dx,dy,dz] = dst_leaf.mesh(RES_PER_NODE, INTERP_TYPE);
                tmpval = qdata.interp_points(src_tree, xx, yy, zz, INTERP_TYPE);
                qdata.set_node_val(dst_leaf, tmpval, RES_PER_NODE);
            end
        end

        %/* ************************************************** */
        function valq = interp_points(src_tree,xq,yq,zq,INTERP_TYPE)
            global RES_PER_NODE;
            src_leaves  = src_tree.leaves();
            valq = zeros(size(xq));
            for src_lvcnt =1:length(src_leaves)
                src_leaf = src_leaves{src_lvcnt};
                indices = qdata.points_in_node(src_leaf, xq, yq);
                if ~any(any(indices)), continue; end;
                xx = xq(indices);
                yy = yq(indices);
                if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                    global CHEB_IMPL
                    if strcmp(CHEB_IMPL, 'IAS')
                        vv = zeros(size(xx));
                        w = src_leaf.data.values;
                        [xmin,xmax,ymin,ymax] = src_leaf.corners;
                        xs = (xx - xmin)*2/(xmax-xmin)-1.0;
                        ys = (yy - ymin)*2/(ymax-ymin)-1.0;
                        for xindx =1:size(xx,1)
                            vv(xindx) =  cheb.chebeval2(w, xs(xindx), ys(xindx));
                        end
                    elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                        w = src_leaf.data.values;
                        vv = w(xx,yy);
                    end
                else
                    [xxr,yyr,zzr,dx,dy,dz] = src_leaf.mesh(RES_PER_NODE);
                    interp_data = src_leaf.data.values;
                    vv = interp2(xxr,yyr,interp_data,xx,yy, INTERP_TYPE);
                end
                valq(indices) = vv;
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
                valt = qdata.interp_points(leaf,xxcc,yycc,zzcc,INTERP_TYPE);

                diff = vale - valt;
                [err, indx] = max(abs(diff(:)));
                [i_row, i_col] = ind2sub(size(diff),indx);
                if err > max_err
                    vale;
                    valt;
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
            qdata.set_node_val(node, fnval, resPerNode);
        end

        %/* ************************************************** */
        function set_node_val(node, fnval, resPerNode)
            global INTERP_TYPE
            if strcmp(INTERP_TYPE, 'CHEBYSHEV')
                w = qdata.get_node_cheb_interpolant(node, fnval, resPerNode);
                node.data.values = w;
            else
                node.data.values = fnval;
            end
            node.data.dim = 1;
            node.data.resolution = resPerNode;
        end

        %/* ************************************************** */
        function [w] = get_node_cheb_interpolant(node, fn_val, resPerNode)
            global CHEB_IMPL;
            if strcmp(CHEB_IMPL, 'IAS')
                w = cheb.chebcoeff(fn_val);
            elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                [xmin xmax ymin ymax] = node.corners;
                w = chebfun2(fn_val, [xmin xmax ymin ymax]);
            end
        end

        %/* ************************************************** */
        function [val] = get_node_values_cheb(node, resPerNode)
            global CHEB_IMPL;
            if strcmp(CHEB_IMPL, 'IAS')
                x = cheb.chebnodes1(resPerNode);
                y = x;
                w = node.data.values;
                val = cheb.chebeval2(w,x,y);
            elseif strcmp(CHEB_IMPL, 'CHEBFUN')
                w = node.data.values;
                [xx,yy,zz,dx,dy,dz] = node.mesh(resPerNode, ...
                                                'CHEBYSHEV');
                val = w(xx,yy);
            end
        end

        %/* ************************************************** */
        function [xx,yy,vv] = grid_points(tree)
            global INTERP_TYPE;
            global RES_PER_NODE;
            xx = []; yy = [];
            cleaves = tree.leaves();
            for lvcnt = 1:length(cleaves)
                cleaf = cleaves{lvcnt};
                if isempty(cleaf.data), continue; end;
                [xr,yr,zr,dx,dy,dz] = cleaf.mesh(RES_PER_NODE,INTERP_TYPE);
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
                    vals = qdata.get_node_values_cheb(cleaf, resPerNode);
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
            [txx,tyy] = qdata.grid_points(tree);
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
            [txx,tyy]   = qdata.grid_points(tree);
            [tvv]       = qdata.grid_data(tree);
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
