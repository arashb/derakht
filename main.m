function main()

clear all;
clear globals; % constants  and preamble

addpath src
addpath src/semilag
addpath src/tree
addpath src/common

% GLOBAL VARIABLES
global verbose;
global gvfreq;
global dim;
global om;
global DEBUG;
global resPerNode;
global maxErrorPerNode;
global maxLevel;
global INTERP_TYPE;

VPREVTSTEP  = 1;
VCURTSTEP   = 2;
VNEXTSTEP   = 3;
V2NEXTSTEP  = 4;

% MAIN SCRIPT
fconc_exact   = @get_gaussian;
fvel_exact    = @get_vel_exact;
fvelx_exact   = @get_velx_exact;
fvely_exact   = @get_vely_exact;
fdo_refine    = @do_refine;         % refinement criterion

% SIMULATION PARAMETERS BASED ON INITIAL CCONC. TREE
max_level_list = [2];% 3];% 4 5 6];
tn      = 50;
for lvl =1:size(max_level_list,2)
    
    % SIMULATION PARAMETERS
    maxErrorPerNode = 1e-30;                             % Error per box
    maxLevel        = max_level_list(lvl);               % Maximum tree depth
    resPerNode      = 7;                                % Resolution per Node
    verbose         = false;
    gvfreq          = 0;
    om              = 1;
    dim             = 2;
    DEBUG           = false;
    INTERP_TYPE     = 'cubic';
    
    % TEMPORAL RESOLUTION
    T_INIT   = 0;
    T_FINAL  = 2*pi;
    tn       = tn*2;
    dt       = T_FINAL/tn;
    t        = T_INIT + [-dt 0 dt 2*dt];
    
    fprintf('--> dt: %d\n',dt);
    
    % OUTPUT FOLDER
    timestamp = datestr(now,30);
    dir_name = ['./results/res_',timestamp,'/'];
    mkdir(dir_name);
    
    % CONSTRUCT AND INIT THE TREES (VELOCITY AND CONCENTRAITION FEILDS)
    fprintf('-> init the trees\n');
    
    % VELOCITY (TIME-INDEPENDENT)
    u = qtree;
    u.insert_function(fvelx_exact,fdo_refine,t(2));
    tree_data.init_data(u,fvelx_exact,resPerNode,t(2));
    
    v = qtree;
    v.insert_function(fvely_exact,fdo_refine,t(2));
    tree_data.init_data(v,fvely_exact,resPerNode,t(2));
    
    % CONCENTRATION TREE
    cinit = qtree;
    cinit.insert_function(fconc_exact,fdo_refine,T_INIT);
    tree_data.init_data(cinit,fconc_exact,resPerNode,T_INIT);
    
    plot_tree(cinit,0);
    
    % COMPUTE INITIAL TREES' ERRORS
    err = zeros(tn+1,2);
    err(1,1)= 0;
    err(1,2)= tree_data.compute_error(cinit, fconc_exact, T_INIT, INTERP_TYPE);
    
    c = cinit;
    for tstep =1:tn
        fprintf('======================================\n');
        
        % ONE STEP SL ADVECTION
        tic;
        cnext = advect_tree_semilag(c,u,v,t,fdo_refine,fconc_exact,fvel_exact);
        toc
        
        % PLOT THE RESULT
        plot_tree(cnext,tstep);
        
        % COMPUTE THE ERROR
        e = tree_data.compute_error(cnext, fconc_exact, t(VNEXTSTEP), INTERP_TYPE);
        err(tstep+1,1) = tstep;
        err(tstep+1,2) = e;
        format longE
        disp(['TN: ', num2str(tstep), '   Error: ', num2str(e,'%12.2e')]);
        
        % PREPARE FOR THE NEXT STEP
        c = cnext;
        t = t + dt;
    end % for time step
    
    % SAVE THE ERROR
    fileID = fopen([dir_name,'error.txt'],'w');
    fprintf(fileID,'%6s %12s\n','tn','LINF');
    fprintf(fileID,'%6d %12.2e\n',err');
    fclose(fileID);
end

% PLOT THE RESULTS
    function plot_tree(cnext, tstep)
        f = figure('Name','SEMI-LAG ADVECTION','Visible','off');
        cnext.plottree(0.5);
        tree_data.plot_data(cnext);
        colorbar;
        
        title(['tstep: ',num2str(tstep)]);
        b=sprintf('%4.3d',tstep);
        s_fig_name = [dir_name,'fig_',b]; % ,
        saveas(f,s_fig_name,'png');
        close(f)
    end

%/* ************************************************** */
    function [val, funcval] = do_refine(qtree,func,t)
        [val, funcval] = tree_do_refine(qtree,func,maxErrorPerNode,maxLevel,resPerNode,t);
    end

end

%/* ************************************************** */
function value = get_gaussian(t,x,y,z)
global om;
xc = 0.5;
yc = 0.5;
xci = 0.6;
yci = 0.5;

[alpha,RHO] = cart2pol(xci-xc,yci-xc);
alphat = alpha + t*om;

[xct,yct] = pol2cart(alphat,RHO);
xct = xct + xc;
yct = yct + yc;

theta = 0;
sigmax = 0.06;
sigmay = 0.06;
value = gaussian(x,y,xct,yct,theta,sigmax,sigmay);
end

%/* ************************************************** */
function value = get_gaussian_2(t,x,y,z)
global om;
xc = 0.5;
yc = 0.5;
theta = 0;
sigmax = 0.05;
sigmay = 0.05;

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
function value = get_zalesak(t,x,y,z)
value = zalesak(x, y);
end

%/* ************************************************** */
function [u,v,w] = get_vel_exact(t,x,y,z)
xc = 0.5;
yc = 0.5;
zc = 0.5;
[u,v,w] = vel_rot(t,x,y,z,xc,yc,zc);
end

%/* ************************************************** */
function u = get_velx_exact(t,x,y,z)
[u,v,w] = get_vel_exact(t,x,y,z);
end

%/* ************************************************** */
function v = get_vely_exact(t,x,y,z)
[u,v,w] = get_vel_exact(t,x,y,z);
end
