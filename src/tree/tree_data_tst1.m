function tree_data_tst1()
% test the creation of quadtree based on inserted function(s)
clear; clear globals; dim=2;  % constants  and preamble
addpath('../common/');
global verbose;
global gvfreq;

% RUN PARAMETERS
maxErrorPerNode = 0.01;      % Error per box
maxLevel        = 20;          % maximum tree depth
verbose         = false;
resPerNode      = 10;
gvfreq          = 1;

% MAIN SCRIPT

% create the first tree
o = qtree;
o.insert_function(@func1,maxErrorPerNode,maxLevel,resPerNode);

% init the fist tree's data with a given function
tree_data.init_data(o,@func3,resPerNode)

% create the second tree 
q = qtree;
q.insert_function(@func2,maxErrorPerNode,maxLevel,resPerNode)

% interpolate second tree's data from the first tree
tree_data.interp(o,q);

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

subplot(3,2,2)
tree_data.plot_grid(q)

subplot(3,2,3)
tree_data.plot_values(o,1)

subplot(3,2,4)
tree_data.plot_values(q,1)

subplot(3,2,5)
tree_data.plot_values(o,2)

subplot(3,2,6)
tree_data.plot_values(q,2)

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

    function value = func3(x,y)
        t = 0;
        z = 0;
        value = zeros(size(x));
        xc = 0.5*ones(size(x)); 
        yc = xc; zc = xc;
        [u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
        value(:,:,:,1) = u;
        value(:,:,:,2) = v;
    end
end

