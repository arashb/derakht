function tst1()
%tst5 Test the interpolation of values between trees with differnt depths
clear; clear globals; % constants  and preamble
addpath tree
addpath semilag
addpath common

global verbose;

% RUN PARAMETERS
maxErrorPerNode = 0.001;                                % Error per box
maxLevel        = 20;                                   % Maximum tree depth
resPerNode      = 15;                                   % Resolution per Node
verbose         = true;

interpMaxErrorPerNode = maxErrorPerNode*0.01;           % Interplant Error per box

% MAIN SCRIPT
fconc = @func1;

c = qtree;
c.insert_function(fconc,interpMaxErrorPerNode,maxLevel,resPerNode);
tree_init_data(c, fconc, resPerNode);

cinterp = qtree;
cinterp.insert_function(fconc,maxErrorPerNode,maxLevel,resPerNode);

% INTERPOLATION
tree_interp(c, cinterp);

% PLOTTING
MS='MarkerSize';

subplot(2,2,1)
c.plottree;
hold on;
[cxx,cyy,cvv] = tree_griddata(c);
plot(cxx(:),cyy(:),'ro',MS,1); hold on;
axis off; axis equal;
title('c');
 
subplot(2,2,2)
cinterp.plottree;
hold on;
[ctxx,ctyy,ctvv] = tree_griddata(cinterp);
plot(ctxx(:),ctyy(:),'ro',MS,1); hold on;
axis off; axis equal;
title('interp. c');

subplot(2,2,3);
plot3(cxx,cyy,cvv,'.',MS,1);
axis off; axis equal;
title('c');

subplot(2,2,4);
%scatter(ctxx,ctyy,ctvv,'.','MarkerFaceColor',ctvv);
plot3(ctxx,ctyy,ctvv,'.',MS,1);
axis off; axis equal;
title('interp. c');
end

%/* ************************************************** */
function value = func1(x,y)
xc = 0.5;
yc = 0.5;
theta = pi*0.5;
sigmax = 0.05;
sigmay = 0.09;
value = gaussian(x,y,xc,yc,theta, sigmax, sigmay);
end
