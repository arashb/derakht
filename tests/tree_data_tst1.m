function tree_data_tst1()
% test the interpolation of multidimensional values between two quadtrees with arbitrary
% tree depth
    clear all;
    clear globals;  % constants  and preamble

    addpath('../src/common/');
    addpath('../src/tree/');

    global verbose;
    global gvfreq;
    global dim;
    global om;
    global resPerNode;
    global maxErrorPerNode;
    global maxLevel;
    global INTERP_TYPE;
    global CHEB_IMPL;
    %global CHEB_KIND;

    % RUN PARAMETERS
    maxErrorPerNode = 1e-4;       % Error per box
    maxLevel        = 15;           % maximum tree depth
    resPerNode      = 3;
    verbose         = false;
    gvfreq          = 0;
    om              = 1;
    dim             = 2;

    %INTERP_TYPE     = 'cubic';
    INTERP_TYPE     = 'CHEBYSHEV';
    CHEB_IMPL       = 'CHEBFUN';
    CHEB_IMPL       = 'IAS';

    % MAIN SCRIPT
    % create the first tree
    disp('first tree')
    o = qtree;
    o.insert_function(@func1, @do_refine);
    qdata.init_data(o, @func1, resPerNode);

    figure
    o.plottree;
    qdata.plot_data(o,1);
    title('initial first value');

    max_err = qdata.compute_error(o, @func1, 0, INTERP_TYPE)


    q = qtree;
    q.insert_function(@func1, @do_refine);
    qdata.interp_tree(o, q);

    figure
    q.plottree;
    qdata.plot_data(q,1)
    title('interpolant first value');

    max_err = qdata.compute_error(q, @func1, 0, INTERP_TYPE)

    % p = qtree;
    % p.insert_function(@func1, @do_refine);
    % qdata.interp_tree(q, p);

    % figure
    % p.plottree;
    % qdata.plot_data(p,1)
    % title('interpolant first value');
    % max_err = qdata.compute_error(p, @func1, 0, INTERP_TYPE)

    % depth=find_depth(o);
    % fprintf('tree depth is %d\n', depth);

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
