clear all; clear globals; dim=2;  % constants  and preamble
global verbose;

% RUN PARAMETERS
maxErrorPerNode = 0.0001;      % Error per box
maxLevel        = 20;          % maximum tree depth
verbose         = false;


% MAIN SCRIPT
func = @gaussian;

% plot the function
res = 1000;
x = linspace(0,1,res);
y = linspace(0,1,res);
[xx,yy] = meshgrid(x,y);
Z = func(xx,yy); 
contourf(xx,yy,Z);
colorbar;
axis off;
hold on;

% create and plot the tree
o = qtree;
o.insert_function(func,maxErrorPerNode,maxLevel);
o.plottree;
axis off; hold on;

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


