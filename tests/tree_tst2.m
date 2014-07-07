function tree_tst2()
% test the reconstruction of the quadtree based on the leaves morton ids
clear; clear globals; dim=2;  % constants  and preamble
addpath('../src/common/');
addpath('../src/tree/');
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
Z1 = func1(0,xx,yy,0);
Z2 = func2(0,xx,yy,0);
Z3 = func3(0,xx,yy,0);

a = qtree;
a.insert_function(@func1,@do_refine);
a.insert_function(@func2,@do_refine);
a.insert_function(@func3,@do_refine);

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

    function value = func1(t,x,y,z)
        xc = 0.75;
        yc = 0.75;
        value = gaussian(x,y,xc,yc);
    end

    function value = func2(t,x,y,z)        
        xc = 0.25;
        yc = 0.25;
        value = gaussian(x,y,xc,yc);
    end

    function value = func3(t,x,y,z)        
        value = gaussian(x,y);
    end

    function val = do_refine(qtree,func,t)
        val = tree_do_refine(qtree, func, maxErrorPerNode, maxLevel, resPerNode,t);
    end
end