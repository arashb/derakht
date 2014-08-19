function [] = rotating_cone_tst1()
clear; close all;
clear global;
addpath('../src/common/');
addpath('../src/semilag/');


VF_TYPE  = 1; % type of velocity field
CF_TYPE  = 1; % type of intial concentration field
ERR_TYPE = 1; % type of error computation (L_2, L_inifinite, ...)

PLOT_SOL   = true;
RES_PATH   = './'; % path to save results
fig_format = '.pdf';

%INTERP_LIST = {'linear','cubic','spline'};
%INTERP_LIST = {'linear'};
%INTERP_LIST  = {'cubic'};
INTERP_LIST = {'spline'};

global gvfreq;
global dim;
global verbose;
global om;
global DEBUG;

dim     = 2;
verbose = false;
DEBUG   = true;
om      = 1.0;

xi      = 0;
xf      = 1;
n_level = 4:8;
n_list  = 2.^n_level;

ti      = 0;
tf      = 2*pi;

fcnt    = 1;
gvfreq  = 0;


for interptypecnt=1:length(INTERP_LIST)
    INTERP_TYPE = INTERP_LIST{interptypecnt};

    dt      = pi/20;
    tn      = (tf - ti)/dt;

for ncnt =1:length(n_list)      % spatial iteration
    n  = n_list(ncnt);
    dx = (xf - xi)/n_list(ncnt);

    fprintf('*************************************\n');
    fprintf('interpolation type: %s\n',INTERP_TYPE);
    fprintf('ti: %d\ntf: %d\ndt: %d\ntn:%d\n',ti,tf,dt,tn);
    fprintf('n: %d\ndx: %d\n',n,dx);

    x = linspace(xi,xf,n+1);
    [xx, yy, zz] = meshgrid(x, x, 1:1);
    height = 1; %round(2*(n+1)/3);      % used for slicing the cylinder

    t = [ti-dt, ti, ti+dt, ti+2*dt];
    [ u, v, w, cinit ] = init_fields(xi, xf, dx, xx, yy, zz, t, VF_TYPE, CF_TYPE);
    csol = cinit;

    crk2 = compute_numerical(cinit, xx, yy, zz, u, v, w, t, dt, tn, INTERP_TYPE, 'rk2');
    error_rk2(interptypecnt,ncnt,fcnt) = compute_error(crk2(3:end-2,3:end-2,:,end), ...
                                                       csol(3:end-2,3:end-2,:), ERR_TYPE);
    fprintf('rk2 error =  %e\n',error_rk2(interptypecnt,ncnt,fcnt));

    c2tl = compute_numerical(cinit, xx, yy, zz, u, v, w, ti, dt, tn, INTERP_TYPE, '2tl');
    error_2tl(interptypecnt,ncnt,fcnt) = compute_error(c2tl(3:end-2,3:end-2,:,end), ...
                                                       csol(3:end-2,3:end-2,:), ERR_TYPE);
    fprintf('2tl error =  %e\n',error_2tl(interptypecnt,ncnt,fcnt));


    if DEBUG
        figure;
       for tscnt =1:tn
        subplot(3,2,1);
        surf(cinit(:,:,height)); title('INITIAL');
        %az = 0; el = 90; view(az, el);
        colorbar;
        %axis off; axis equal; drawnow;
        
        subplot(3,2,2);
        surf(csol(:,:,height)); title('ANALYTICAL');
        %az = 0; el = 90; view(az, el);
        colorbar;
        %axis off; axis equal; drawnow;
        
        subplot(3,2,3)
        surf(crk2(:,:,height,tscnt+1)); title(['RK2 ts:',num2str(tscnt)]);
        %az = 0; el = 90; view(az, el);
        colorbar;
        %axis off; axis equal; drawnow;

        subplot(3,2,4);
        surf(c2tl(:,:,height,tscnt+1)); title(['2TL ts:',num2str(tscnt)]);
        %az = 0; el = 90; view(az, el);
        colorbar;
        %axis off; axis equal; drawnow;

        subplot(3,2,5)
        diff_rk2 = crk2(:,:,:,tscnt+1) - csol(:,:,:);
        surf(diff_rk2(:,:,height)); title('RK2 Error');
        az = 0; el = 90; view(az, el);
        colorbar;
        %axis off; axis equal; drawnow;

        subplot(3,2,6)
        diff_c2l = c2tl(:,:,:,tscnt+1) - csol(:,:,:);
        surf(diff_c2l(:,:,height)); title('2TL Error');
        az = 0; el = 90; view(az, el);
        colorbar;
        %axis off; axis equal; drawnow;
        pause;
        end
    end

    % set temporal resolution manually
    dt = dt*0.5;
    tn = tn*2;

end                             % spational resolution loop 
end                             % interpolation loop

% plot the error
m = {'h','o','*','.','x','s','d','^','v','>','<','p','h'};
ls = {'-','--', ':', '-.'};
cmap = {'r','g', 'b', 'c', 'm', 'y', 'k' , 'w'};
cmapcnt = 1;
legendInfo = {};
s_fig_name = ['error_cmp_', datestr(now)];
f = figure('Name',s_fig_name,'units','normalized','outerposition',[0 0 1 1]);
for interptypecnt=1:length(INTERP_LIST)
        % create and save the error graph
        figure(f)
        h = loglog(1./n_list,error_rk2(interptypecnt,:,fcnt));
        grid on;
        mi = mod(cmapcnt,length(m)) + 1;
        set(h, 'Marker', m{mi});
        lsi = mod(cmapcnt,length(ls)) + 1;
        set(h, 'LineStyle', ls{lsi});
        ci = mod(cmapcnt, length(cmap)) + 1;
        set(h, 'Color', cmap{ci});
        cmapcnt= cmapcnt+1;
        mes = sprintf('rk2-%s',INTERP_LIST{interptypecnt});
        legendInfo{end+1}= mes;
        hold on;
        i = loglog(1./n_list,error_2tl(interptypecnt,:,fcnt));
        grid on;
        mi = mod(cmapcnt,length(m)) + 1;
        set(i, 'Marker', m{mi});
        lsi = mod(cmapcnt,length(ls)) + 1;
        set(i, 'LineStyle', ls{lsi});
        ci = mod(cmapcnt, length(cmap)) + 1;
        set(i, 'Color', cmap{ci});
        cmapcnt= cmapcnt+1;
        mes = sprintf('twotl-%s',INTERP_LIST{interptypecnt});
        legendInfo{end+1}= mes;
        hold on;
end
legend(legendInfo);
title('Error Comparision');
xlabel('Domain/Time Resolution');
ylabel('Error (logarithmic)');
s_fig_path = [RES_PATH,s_fig_name,fig_format];
saveas(f,s_fig_path);
end
