function [] = zalesak_tst1()
% ZALESAK_TST1 check the semi-Lagrangian interpolation accurcy for the Zalesak test case
clear; close all;
clear global;

addpath('../src/common/');
addpath('../src/semilag/');

VF_TYPE     = 1;                    % type of velocity field
CF_TYPE     = 0;                    % type of intial concentration field: 0:Zalesak cylinder 1:Gaussian  
ERR_TYPE    = 1;                    % type of error computation (L_2, L_inifinite, ...)

PLOT_SOL    = 1;
RES_PATH    = './';                % path to save results
fig_format  = '.pdf';

%INTERP_LIST = {'linear','cubic','spline'};
INTERP_LIST = {'spline'};


global gvfreq;
global om;
global dim;
global verbose;

gvfreq  = 0;
dim     = 2;
verbose = false;
om      = 1.0;

n_level = [4];
n_list  = 2.^n_level;
cfl     = 1.0;

% spatial domain
xi      = 0.0;
xf      = 1.0;
dx      = (xf - xi)/n_list(1)

% temporal domain
ti      = 0.0;
tf      = 2*pi;
dt      = cfl*dx/om;
tn      = tf/dt;
%tn      = 10;
%tf      = tn*dt;

fcnt = 1;
ncnt = 1;

n = n_list(1);

fprintf('*************** Simulation Data **********************\n');
fprintf('ti: %d\n tf: %d\n',ti,tf);
fprintf('n: %d\ndx: %d\ndt: %d\n',n,dx,dt);
fprintf('*************************************\n');

x = linspace(xi,xf,n+1);
[xx, yy, zz] = meshgrid(x, x, 1:1);
height = 1;      % used for slicing the cylinder

t = [ti-dt, ti, ti+dt, ti+2*dt];
[ u, v, w, cinit ] = init_fields(xi, xf, dx, xx, yy, zz, t, VF_TYPE, CF_TYPE);

color_spec = {'r','g','b','c'};
legendInfo{1} = 'initial';

f = figure('Name','Zalesak Disk','units','normalized','outerposition',[0 0 1 1]);
contour(cinit(:,:,height),1,'k');
hold on
axis off;
%axis equal;

for interptypecnt=1:length(INTERP_LIST)
    INTERP_TYPE = INTERP_LIST{interptypecnt};
    fprintf('*************************************\n');
    fprintf('interpolation type: %s\n',INTERP_TYPE);

    crk2 = compute_numerical(cinit, xx, yy, zz, u, v, w, t, dt, tn, INTERP_TYPE, 'rk2');
    %c2tl = compute_numerical(cinit, xx, yy, zz, u, v, w, t, dt, tn, INTERP_TYPE, '2tl');

    contour(crk2(:,:,height,end),1,color_spec{interptypecnt});
    legendInfo{end+1} = [['rk2-'],INTERP_TYPE];
end
legend(legendInfo);
%shg
fig_format = '.pdf';
fig_path   = './';
s_fig_name = ['zalesak_disk_cmpr_',datestr(now)];
s_fig_path = [fig_path,s_fig_name,fig_format];
saveas(f,s_fig_path);

end
