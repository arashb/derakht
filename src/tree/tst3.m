function tst3()
% test the creation and merge of several quadtrees
clear all; clear globals; dim=2;  % constants  and preamble
global verbose;

% RUN PARAMETERS
maxErrorPerNode = 0.01;      % Error per box
maxLevel        = 20;          % maximum tree depth
resPerNode      = 10;          % Resolution per box
verbose         = false;

% plot the function
res = 100;
x = linspace(0,1,res);
y = linspace(0,1,res);
[xx,yy] = meshgrid(x,y);
Z1 = func1(xx,yy);
Z2 = func2(xx,yy);
Z3 = gaussian(xx,yy);

% create and plot the tree for func1
a = qtree;
b = qtree;
c = qtree;
a.insert_function(@func1,maxErrorPerNode,maxLevel,resPerNode);
b.insert_function(@func2,maxErrorPerNode,maxLevel,resPerNode);
c.insert_function(@gaussian,maxErrorPerNode,maxLevel,resPerNode);
d = qtree.merge(a,b);
e = qtree.merge(d,c);

subplot(2,2,1);
contour(xx,yy,Z1);
a.plottree;
subplot(2,2,2);
contour(xx,yy,Z2);
b.plottree;
subplot(2,2,3);
contour(xx,yy,Z3);
c.plottree;
subplot(2,2,4);
e.plottree;
axis off; hold on;

if verbose
    % print morton ids, all nodes
    %disp(' all nodes');
    %c.print_mids;
    % print morton ids, leaves only
    disp('  leaves only');
    a.print_mids(true);
end
depth=find_depth(a);
fprintf('tree depth is %d\n', depth);


    function value = func1(x,y)
        xc = 0.75;
        yc = 0.75;
        value = gaussian(x,y,xc,yc, 0, 0.02, 0.05);
    end

    function value = func2(x,y)        
        xc = 0.25;
        yc = 0.25;
        value = gaussian(x,y,xc,yc);
    end

end
