function tree_data_tst1()
% test the interpolation of multidimensional values between two quadtrees with arbitrary
% tree depth
    clear;
    clear globals;  % constants  and preamble

    addpath('../src/common/');
    addpath('../src/tree/');
    % addpath('../src/cheb/');

    global verbose;
    global gvfreq;
    global dim;
    global om;
    global resPerNode;
    global maxErrorPerNode;
    global maxLevel;
    global INTERP_TYPE;

    % RUN PARAMETERS
    maxErrorPerNode = 0.0001;      % Error per box
    maxLevel        = 20;          % maximum tree depth
    resPerNode      = 14;
    verbose         = false;
    gvfreq          = 0;
    om              = 1;
    dim             = 2;
    % INTERP_TYPE     = 'cubic';
    INTERP_TYPE     = 'CHEBYSHEV';

    % MAIN SCRIPT
    % create the first tree
    o = qtree;
    o.insert_function(@func1,@do_refine);
    tree_data.init_data(o,@func1,resPerNode)

    % TODO: CHECK THE ERROR

    % create the second tree
    q = qtree;
    q.insert_function(@func2,@do_refine)

    % interpolate second tree's data from the first tree
    % tree_data.interp(o,q);

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

    % PLOTTING
    subplot(3,2,1)
    tree_data.plot_grid(o)
    title('initial tree');

    subplot(3,2,2)
    tree_data.plot_grid(q)
    title('interpolant tree');

    subplot(3,2,3)
    o.plottree;
    tree_data.plot_data(o,1)
    title('initial first value');

    subplot(3,2,4)
    q.plottree;
    tree_data.plot_data(q,1)
    title('interpolant first value');

    % subplot(3,2,5)
    % o.plottree;
    % tree_data.plot_data(o,2)
    % title('initial second value');

    % subplot(3,2,6)
    % q.plottree;
    % tree_data.plot_data(q,2)
    % title('interpolant second value');

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
        value = zeros(size(x));
        xc = 0.5*ones(size(x));
        yc = xc; zc = xc;
        [u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
        value(:,:,:,1) = u;
        value(:,:,:,2) = v;
    end

    function val = do_refine(qtree,func,t)
        val = tree_do_refine(qtree, func, maxErrorPerNode, maxLevel, resPerNode,t);
    end
end
