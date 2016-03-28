function plot_mt()
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
    maxErrorPerNode = 1e-4;      % Error per box
    maxLevel        = 20;          % maximum tree depth
    resPerNode      = 8;          % Resolution per box
    verbose         = false;
    INTERP_TYPE     = 'CUBIC';
    dim=2;

    % MAIN SCRIPT
    f1 = @func1;
    f2 = @func2;
    f3 = @func3;

    np = 4;

    % create and plot the trees
    c = qtree;
    c.insert_function(f2,@do_refine);
    c.plottree_tikzNP('c.tex',np)

    v = qtree;
    v.insert_function(f1,@do_refine);
    v.plottree_tikzNP('v.tex',np)

    %  COMPLETE MERGE
    cm = qtree.merge(c,v);
    cm.plottree_tikzNP('cm.tex',np);

    sl = qtree.semi_merge(c,v,np);
    c.plottree_tikzK('csm.tex',sl);
    v.plottree_tikzK('vsm.tex',sl);


    if verbose
        % print morton ids, all nodes
        disp(' all nodes');
        o.print_mids;
        % print morton ids, leaves only
        disp('  leaves only');
        o.print_mids(true);
    end

    function value = func1(t,x,y,z)
        xc = 0.7;
        yc = 0.7;
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
