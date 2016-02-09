function tree_tst1()
% test the creation of quadtree based on inserted function(s)
    close all; clear all; clear globals;   % constants  and preamble
    addpath('../src/common/');
    addpath('../src/tree/');

    % GLOBAL VARIABLES
    global verbose;
    global gvfreq;
    global dim;
    global om;
    global resPerNode;
    global maxErrorPerNode;
    global maxLevel;
    global INTERP_TYPE;
    global vis;

    % RUN PARAMETERS
    maxErrorPerNode = 0.0001;      % Error per box
    maxLevel        = 20;          % maximum tree depth
    resPerNode      = 10;          % Resolution per box
    verbose         = false;
    INTERP_TYPE     = 'CUBIC';
    dim=2;

    % MAIN SCRIPT
    f1 = @func1;
    f2 = @func2;
    f3 = @func3;

    % plot the function
    res = 1000;
    x = linspace(0,1,res);
    y = linspace(0,1,res);
    [xx,yy] = meshgrid(x,y);

    Z1 = f1(0,xx,yy,0);
    Z2 = f2(0,xx,yy,0);
    Z3 = f3(0,xx,yy,0);

    hold on
    contour(xx,yy,Z1);
    contour(xx,yy,Z2);
    contour(xx,yy,Z3);
    colorbar;
    axis off;


    % create and plot the tree
    o = qtree;
    o.insert_function(f1,@do_refine);
    o.insert_function(f2,@do_refine);
    o.insert_function(f3,@do_refine);
    hold on;
    o.plottree;
    axis off;

    if verbose
        % print morton ids, all nodes
        disp(' all nodes');
        o.print_mids;
        % print morton ids, leaves only
        disp('  leaves only');
        o.print_mids(true);
    end
    depth=find_depth(o);
    fprintf('tree depth is %d\n', depth);

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
        xc = 0.5;
        yc = 0.5;
        value = gaussian(x,y,xc,yc);
    end

    function val = do_refine(qtree,func,t)
        val = tree_do_refine(qtree, func, maxErrorPerNode, maxLevel, resPerNode,t);
    end
end
