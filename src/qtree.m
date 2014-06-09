classdef qtree < handle
% implements a quatree data structure
%    point based construction in 2D.
%    preorder and postorder traversals
%      
% see "doc qtree" for documentation on methods and classes
% see "nbody.m" for the 

  
%/* ************************************************** */
properties
  kids    % kid structure  kids.m (size) kids.ptr{i}  pointeres
  parent  % parent of the node
  level   % level of the node 
  data    % averaging-related data
  isleaf  % true or false
  anchor  % [x;y] coordinates of lower left point of a node
end

methods
%/* ************************************************** */
function this=qtree(parent,level,anchor)
% function this=qtree(parent,level,anchor)
% constructor for qtree 

if nargin<1, parent = [];           end;
if nargin<2, level  = 0;            end;
if nargin<3, anchor = [0,0];        end;

this.kids   = [];
this.parent = parent;
this.level  = level;
this.data   = [];
this.isleaf = true;
this.anchor = anchor;

end


%/* ************************************************** */
function h=width(this)
% function h=width(this)
% get the width of the square that corresponds to this node
  h = 1/2^this.level;
end

%/* ************************************************** */
function [xmin,xmax,ymin,ymax]=corners(this)
% function [xmin,xmax,ymin,ymax]=corners(this)  
%  get the   coordinates of the corner points of a node
% lower left point [xmin; ymin], upper left point [xmin; ymax] etc. 
  
  h = width(this);
  xmin = this.anchor(1);
  ymin = this.anchor(2);

  xmax = xmin + h;
  ymax = ymin + h;

% to capture points on domain boundaries, modify xmax,ymax
  if xmax==1, xmax=1+eps; end 
  if ymax==1, ymax=1+eps; end;
end
  
  
%/* ************************************************** */
function create_kids(this)  
% function create_kids(this)    
% create the kid datastructure of this node

assert(isempty(this.kids),'Kids already created');

kidlevel= this.level+1;
kidh    = this.width/2;

% X Y
% 0 0 child  ( lower left  )
% 0 1 child  ( lower right )
% 1 0 child  ( upper left  )
% 1 1 child  ( upper right )
kids=[ [0 0]; [1 0]; [ 0 1]; [1 1] ];
anchor= zeros(4,2);
for i=1:4
  anchor(i,:) = this.anchor + kidh * kids(i, :);
  this.kids{i} = qtree(this, kidlevel,anchor(i,:));
end

%  update parent information
this.isleaf = false;
end

%/* ************************************************** */
% points = dim x num_of_points
function insert_function(this, func, maxErrPerNode, maxLevel)
% This is the main construction routine for qtree.
 
if this.isleaf
    if ~do_refine(this, func, maxErrPerNode, maxLevel), return; end
    % this node will be split;
    create_kids(this);
end

% now we insert the function to the kids
for k=1:4
    this.kids{k}.insert_function( func, maxErrPerNode, maxLevel);
end

    %/* ************************************************** */
    function [val] = do_refine(this, func, maxErrPerNode, maxLevel)
        %global verbose
        if this.level == maxLevel, val = false; return; end;
        m = 15;
        err = compute_error_rg(this, func, m);
        %err= compute_error_box_center(this, func)
%         if verbose
%             fprintf('node  at level %d: anchor:[%1.4f %1.4f] error: %e\n',...
%                 this.level,this.anchor(1),this.anchor(2), err);
%         end
        if err <=  maxErrPerNode, val = false; return; end;
        val = true;
    end
    
    %/* ************************************************** */
    function err= compute_error_rg(this, func, m)
        % get the coordinates of the corners of the current box
        [xmin,xmax,ymin,ymax]=corners(this);
        
        % create the regular grid in the box with resolution (m+1)^2
        xr = linspace(xmin, xmax, m+1);
        dx = (xmax - xmin)/m;
        yr = linspace(ymin, ymax, m+1);
        dy = (ymax - ymin)/m;
        [xxr, yyr] = meshgrid(xr,yr);
        
        % compute the function values on the regular grid points
        fre = func(xxr,yyr);
        
        % compute the center of the regular grid cells
        xxc = xxr+dx/2;
        yyc = yyr+dy/2;
        xxcc = xxc(1:end-1,1:end-1);
        yycc = yyc(1:end-1,1:end-1);
        
        % compute the exact values on the centers
        fce = func(xxcc, yycc);
        
        % interpolate the function values on the center points
        fci = interp2(xxr, yyr, fre, xxcc, yycc);
        
        % compute interpolation error
        diff = fci - fce;
        err = max(max(abs(diff)));
    end
    
    %/* ************************************************** */
    function err= compute_error_box_center(this, func)
        % get the coordinates of the corners of the current box
        [xmin,xmax,ymin,ymax]=corners(this);
        
        % create the regular grid in the box with resolution (m+1)^2
        xr = linspace(xmin, xmax, 2);
        %dx = (xmax - xmin)/m;
        yr = linspace(ymin, ymax, 2);
        %dy = (ymax - ymin)/m;
        [xxr, yyr] = meshgrid(xr,yr);
        
        % compute the function values on the regular grid points
        fre = func(xxr,yyr);
        
        % compute the center of the regular grid cells
        xc = 0.5*(xmin + xmax);
        yc = 0.5*(ymin + ymax);
        
        % compute the exact values on the centers
        fce = func(xc, yc);
        
        % interpolate the function values on the center point
        fci = interp2(xxr, yyr, fre, xc, yc, 'linear');
        
        % compute interpolation error
        diff = fci - fce;
        err = abs(diff);
    end
    
end

%/* ************************************************** */
function  preorder(this, visit, prune, user_data)
% function  preorder(this, visit, prune, user_data)  
% preorder traversal with prunning  
% visit(node,data)
% prune(node,data) can be empty
% user_data: user data structure for calculations

doPrune = false;
if ~isempty(prune),   doPrune= prune(this, user_data);  end
visit(this,user_data);
if ~this.isleaf & ~doPrune
   for k=1:4
     this.kids{k}.preorder( visit, prune, user_data);
   end
end
end


%/* ************************************************** */
function  postorder(this, visit, prune, user_data)
% function  postorder(this, visit, prune, user_data)  
% preorder traversal with prunning  
% visit(node,data)
% prune(node,data) can be empty
% user_data: user data structure for calculations

doPrune = false;
if ~isempty(prune),   doPrune= prune(this, user_data);  end
if ~this.isleaf & ~doPrune
   for k=1:4
     this.kids{k}.postorder( visit, prune, user_data);
   end
end
visit(this,user_data);
end

%/* ************************************************** */
function list=leaves(this)
% function list=leaves(this)  
% collects all the leaves in a single array  
  list={};
  cnt=0;
  function visit(this,dummy)
    if this.isleaf
      cnt=cnt+1;
      list{cnt}=this;
    end
  end
  this.preorder(@visit,[],[]);
end

%/* ************************************************** */
function list=nodes(this)
% function list=nodes(this)  
% collects all the nodes in a single array  
  list={};
  cnt=0;
  function visit(this,dummy)
      cnt=cnt+1;
      list{cnt}=this;
  end
  this.preorder(@visit,[],[]);
end

%/* ************************************************** */
function print(this)
% function print(this)  
% print information about a node  
  print_node = @(this,dummy)...
      fprintf('node  at level %d: anchor:[%1.4f %1.4f]\n',...
              this.level,this.anchor(1),this.anchor(2));
  this.preorder(print_node,[],[]);
end

%/* ************************************************** */
function print_mids(this,leaves_only)
% function print_mids(this,leaves_only)  
% print the Morton IDs of nodes (using the class morton_id)
% leaves_only : if true, it only prints the leaves.
  mid = morton_id;
  if nargin<2, leaves_only=false; end;
  function print_node(this,dummy)
    if leaves_only & ~this.isleaf, return; end;
    fprintf('mid:');
    id = mid.id(this.level,this.anchor);
    %mid.print(id);
    fprintf(' %20u at level %2d: anchor:[%1.4f %1.4f]\n',...
            id, this.level,this.anchor(1),this.anchor(2));

  end
  this.preorder(@print_node,[],[]);    
end

%/* ************************************************** */
function depth=find_depth(this)
% function depth=find_depth(this)  
% function depth=find_depth(this)  
% depth : depth of the tree.  
   depth = 0;
   function v=fdt(this,dummy)
      depth=max(this.level,depth);
   end
   this.preorder(@fdt,[],[]);
end

%/* ************************************************** */
function  plotnode(this,width) 
% function  plotnode(this,width)   
% plot the rectangul corresponding to the node  
  if nargin<2,width=1;end;
  rectangle('Curvature',[0.00,0.00],...
     'Position', [this.anchor(1),this.anchor(2),1/2^this.level,1/2^this.level], 'LineWidth',width);
end

%/* ************************************************** */
function plottree(this,markersize)
% function plottree(this,markersize)  
%  function plottree(this,markersize)
%  plots the whole tree.  
  hold on;
  if nargin<2,markersize=2;end;
  plotnode=@(node,dummy)node.plotnode(markersize);
  this.preorder(plotnode,[],[]);
  hold off;
end
end % methods

methods (Static)
%/* ************************************************** */
function lv_mids = tree2mids(tree, leaves_only)
    if nargin < 2, leaves_only = false; end;
    if leaves_only
        lv_list = tree.leaves();
    else
        lv_list = tree.nodes();
    end
    lv_mids = zeros(size(lv_list),'uint64');
    mid = morton_id;
    for i=1:length(lv_list)
        qt = lv_list{i};
        lv_mids(i) = mid.id(qt.level,qt.anchor);
    end
end

%/* ************************************************** */
function tree = mids2tree(mids_list, leaves_only)
    if nargin < 2, leaves_only = false; end;
    if leaves_only
    else
        [tree, counter] = construct_tree([], mids_list, 1);
    end
    
    %/* ************************************************** */
    function [node, counter] = construct_tree(parent, mids_list, counter)
        %if counter > length(mids_list), return; end;
        global verbose
        id = mids_list(counter);
        mido = morton_id;
        [lvl, anc] = mido.id2node(id);
        if verbose
            fprintf('construncting: %20u at level %2d: anchor:[%1.4f %1.4f]\n',...
                id, lvl,anc(1),anc(2));
        end
        node = qtree(parent, lvl, anc);
        counter = counter+1;
        
        % check for termination
        if counter > length(mids_list), return; end;
        [nlvl, nanc] = mido.id2node(mids_list(counter));
        if nlvl <= lvl, return; end;
        
        node.isleaf = false;
        % construnct the subtrees
        [node.kids{1}, counter] = construct_tree(node, mids_list, counter);
        [node.kids{2}, counter] = construct_tree(node, mids_list, counter);
        [node.kids{3}, counter] = construct_tree(node, mids_list, counter);
        [node.kids{4}, counter] = construct_tree(node, mids_list, counter);
    end
end

%/* ************************************************** */
function treec = merge(treea,treeb)
    lvs_ids1 = qtree.tree2mids(treea);
    lvs_ids2 = qtree.tree2mids(treeb);
    lvs_ids_merged = unique([lvs_ids1, lvs_ids2]);    
    treec = qtree.mids2tree(lvs_ids_merged);
end
end % methods static


end % classdef
