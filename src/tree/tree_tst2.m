function tree_tst2()
% test the reconstruction of the quadtree based on the leaves morton ids
clear; clear globals; dim=2;  % constants  and preamble
global verbose;

% RUN PARAMETERS
maxErrorPerNode = 0.001;      % Error per box
maxLevel        = 20;          % maximum tree depth
resPerNode      = 10;          % Resolution per box
verbose         = true;

res = 100;
x = linspace(0,1,res);
y = linspace(0,1,res);
[xx,yy] = meshgrid(x,y);
Z1 = func1(xx,yy);
Z2 = func2(xx,yy);
Z3 = gaussian(xx,yy);

a = qtree;
a.insert_function(@func1,maxErrorPerNode,maxLevel,resPerNode);
a.insert_function(@func2,maxErrorPerNode,maxLevel,resPerNode);
a.insert_function(@gaussian,maxErrorPerNode,maxLevel,resPerNode);

subplot(1,2,1);
contour(xx,yy,Z1);
hold on;
contour(xx,yy,Z2);
contour(xx,yy,Z3);
a.plottree;
axis off;

if verbose
    % print morton ids, all nodes
    disp(' all nodes');
    a.print_mids();
    % print morton ids, leaves only
%     disp('  leaves only');
%     a.print_mids(true);
end
depth=find_depth(a);
fprintf('tree depth is %d\n', depth);


% get the ids
lvs_ids = qtree.tree2mids(a);
mid = morton_id;
for i =1:length(lvs_ids)
    [lvl, anchor] = mid.id2node(lvs_ids(i));
    if verbose
        fprintf('mid: %20u at level %2d: anchor:[%1.4f %1.4f]\n',lvs_ids(i), ...
        lvl, anchor(1), anchor(2));
    end
end

% get the tree from the ids
aclone = qtree.mids2tree(lvs_ids);
subplot(1,2,2);
contour(xx,yy,Z1);
hold on;
contour(xx,yy,Z2);
contour(xx,yy,Z3);
aclone.plottree;
axis off;

    function value = func1(x,y)
        xc = 0.75;
        yc = 0.75;
        value = gaussian(x,y,xc,yc);
    end

    function value = func2(x,y)        
        xc = 0.25;
        yc = 0.25;
        value = gaussian(x,y,xc,yc);
    end
end