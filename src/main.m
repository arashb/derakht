function main()
clear; clear globals; % constants  and preamble
addpath semilag
addpath tree
addpath common

global verbose;
global gvfreq;
global dim;
global DEBUG;
global resPerNode;
global maxErrorPerNode;
global maxLevel;
global INTERP_TYPE;

% RUN PARAMETERS
maxErrorPerNode = 0.0001;        % Error per box
maxLevel        = 20;           % Maximum tree depth
resPerNode      = 15;          % Resolution per Node
verbose         = true;
gvfreq          = 1;
dim             = 2;
DEBUG           = false;
INTERP_TYPE     = 'cubic';

writerObj = VideoWriter('movie.avi');
open(writerObj);

% MAIN SCRIPT
fconc_exact   = @conc_exact_2;
fvelx   = @velx_exact;
fvely   = @vely_exact;

% CONSTRUCT AND INIT THE TREES VALUES
fprintf('-> init the trees\n');

% CONCENTRATION
cinit = qtree;
tinit = 0;
cinit.insert_function(fconc_exact,@do_refine,tinit);
tree_data.init_data(cinit,fconc_exact,resPerNode,tinit);

% SIMULATION PARAMETERS BASED ON INITIAL CCONC. TREE
cdepth      = cinit.find_depth();
cwidth      = 1/2^cdepth;
dx          = cwidth/resPerNode;
cfl         = 100;
om          = 1;
dt          = cfl*dx/om;
VPREVTSTEP  = 1;
VCURTSTEP   = 2;
VNEXTSTEP   = 3;
V2NEXTSTEP  = 4;
t           = tinit + [-dt 0 dt 2*dt];
tn          = 3;

fprintf('--> init tree depth: %d\n',cdepth);
fprintf('--> cfl: %d\n',cfl);
fprintf('--> dt: %d\n',dt);

c = cinit;
for tstep =1:tn
    fprintf('======================================\n');
    fprintf('time step: %d time: %f \n',tstep, t(VCURTSTEP));
    fprintf('======================================\n');

    cnext = advect_tree_semilag(c,fvelx,fvely,t,@do_refine,fconc_exact,@vel_exact);
    
    plot_res(c,cnext);
    
    c = cnext;
    t = t + dt;    
end

close(writerObj);

    % PLOT THE RESULTS
    function plot_res(c, cnext)        
        f = figure('Name','SEMI-LAG ADVECTION');
        subplot(1,2,1);
        c.plottree(0.5);
        tree_data.plot_data(c)
        colorbar;
        title(['tstep: ',num2str(tstep-1)]);
        
        subplot(1,2,2);
        cnext.plottree(0.5);
        tree_data.plot_data(cnext);
        colorbar;
        
        title(['tstep: ',num2str(tstep)]);
        
%         frame = getframe;
%         writeVideo(writerObj,frame);
        
        s_fig_name = ['results_',num2str(tstep)]; % ,datestr(now)
        saveas(f,s_fig_name,'pdf');
    end

%/* ************************************************** */
    function [val, funcval] = do_refine(qtree,func,t)
        [val, funcval] = tree_do_refine(qtree,func,maxErrorPerNode,maxLevel,resPerNode,t);
    end

end

%/* ************************************************** */
function value = conc_exact_1(t,x,y,z)
xc = 0.5;
yc = 0.5;
om = 1;
theta = -om*t;
sigmax = 0.05;
sigmay = 0.12;
value = gaussian(x,y,xc,yc,theta, sigmax, sigmay);
end

%/* ************************************************** */
function value = conc_exact_2(t,x,y,z)
xc = 0.5;
yc = 0.5;
xci = 0.625;
yci = 0.375;
om = 1;

[alpha,RHO] = cart2pol(xci-xc,yci-xc);
alphat = alpha + t*om;

[xct,yct] = pol2cart(alphat,RHO);
xct = xct + xc;
yct = yct + yc;

theta = 0;
sigmax = 0.05;
sigmay = 0.05;
value = gaussian(x,y,xct,yct,theta,sigmax,sigmay);
end

%/* ************************************************** */
function value = conc_exact_22(t,x,y,z)
xc = 0.5;
yc = 0.5;
theta = 0;
sigmax = 0.05;
sigmay = 0.05;
om = 1;

xc1 = 0.75;
yc1 = 0.75;
[alpha1,RHO1] = cart2pol(xc1-xc,yc1-xc);
alphat1 = alpha1 + t*om;
[xct,yct] = pol2cart(alphat1,RHO1);
xct = xct + xc;
yct = yct + yc;
value1 = gaussian(x,y,xct,yct,theta,sigmax,sigmay);

xc2 = 0.25;
yc2 = 0.25;
[alpha2,RHO2] = cart2pol(xc2-xc,yc2-xc);
alphat2 = alpha2 + t*om;
[xct,yct] = pol2cart(alphat2,RHO2);
xct = xct + xc;
yct = yct + yc;
value2 = gaussian(x,y,xct,yct,theta,sigmax,sigmay);

value = value1 + value2;

end

%/* ************************************************** */
function value = conc_exact_3(t,x,y,z)
xi = 0;
xf = 1;
value = slotted_cylinder( xi, xf, x, y, z);
end

%/* ************************************************** */
function [u,v,w] = vel_exact(t,x,y,z)
    xc = 0.5;
    yc = 0.5;
    zc = 0.5;
    [u,v,w] = vel_rot(0,x,y,z,xc,yc,zc);
end

%/* ************************************************** */
function u = velx_exact(t,x,y,z)
    [u,v,w] = vel_exact(t,x,y,z);
end

%/* ************************************************** */
function v = vely_exact(t,x,y,z)
    [u,v,w] = vel_exact(t,x,y,z);
end




