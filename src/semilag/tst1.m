function [] = tst1()
clear; close all;
clear global;
addpath('../common/');
PLOT_SOL = 1;
VF_TYPE = 1;                    % type of velocity field
CF_TYPE = 1;                    % type of intial concentration field
ERR_TYPE = 1;                   % type of error computation (L_2, L_inifinite, ...)
RES_PATH = './';                % path to save results
fig_format = '.pdf';
INTERP_LIST = {'linear','cubic','spline'};
% INTERP_LIST = {'linear'};
% INTERP_LIST = {'cubic'};
% INTERP_LIST = {'spline'};
VFREQ_LIST = [5];
global gvfreq;
global dim;
global verbose;

dim     = 3;
verbose = false;

xi      = 0;
xf      = 1;
ti      = 0;
n_level = [4 5 6 7];
n_list  = 2.^n_level;
tn      = 1;
cfl     = 1;
om      = 1;                    % this is also defined in vel_rot function
for fcnt = 1:length(VFREQ_LIST)
    gvfreq = VFREQ_LIST(fcnt);
for interptypecnt=1:length(INTERP_LIST)
    INTERP_TYPE = INTERP_LIST{interptypecnt};
    dx = (xf - xi)/n_list(1)/VFREQ_LIST(fcnt);
    dt = cfl/n_list(1)/gvfreq;
for ncnt =1:length(n_list)
    n = n_list(ncnt);
    tf = tn*dt;
    fprintf('*************************************\n');
    fprintf('interpolation type: %s\n',INTERP_TYPE);
    fprintf('ti: %d tf: %d\n',ti,tf);
    fprintf('n: %d\ndx: %d\ndt: %d\nvfreq: %d\n',n,dx,dt,gvfreq);
    x = linspace(xi,xf,n+1);
    [xx, yy, zz] = meshgrid(x, x, x);
    height = round(2*(n+1)/3);      % used for slicing the cylinder
    
    [ u, v, w, cinit ] = init_fields(xi, xf, dx, xx, yy, zz, ti, dt, VF_TYPE, CF_TYPE);    
    csol = compute_analytical(xi, xf, ti, tf, xx, yy, zz, CF_TYPE);
    
    crk2 = compute_numerical(cinit, xx, yy, zz, u, v, w, ti, dt, tn, INTERP_TYPE, 'rk2');
    error_rk2(interptypecnt,ncnt,fcnt) = compute_error(crk2(3:end-2,3:end-2,3:end-2,end), csol(3:end-2,3:end-2,3:end-2), ERR_TYPE);
    fprintf('rk2 error =  %e\n',error_rk2(interptypecnt,ncnt,fcnt));

    c2tl = compute_numerical(cinit, xx, yy, zz, u, v, w, ti, dt, tn, INTERP_TYPE, '2tl');
    error_2tl(interptypecnt,ncnt,fcnt) = compute_error(c2tl(3:end-2,3:end-2,3:end-2,end), csol(3:end-2,3:end-2,3:end-2), ERR_TYPE);
    fprintf('2tl error =  %e\n',error_2tl(interptypecnt,ncnt,fcnt));
  
    
    if PLOT_SOL
        figure;
        subplot(3,2,1);
        surf(cinit(:,:,height)); title('INITIAL'); az = 0; el = 90; view(az, el);
        drawnow;
        colorbar;
        subplot(3,2,2);
        surf(csol(:,:,height)); title('ANALYTICAL'); az = 0; el = 90; view(az, el);
        drawnow;
        colorbar;
        subplot(3,2,3)
        surf(crk2(:,:,height,end)); title('RK2');az = 0; el = 90; view(az, el);
        drawnow
        colorbar;
        subplot(3,2,4);
        surf(c2tl(:,:,height,end)); title('2TL'); az = 0; el = 90; view(az, el);
        drawnow
        colorbar;
        subplot(3,2,5)
        diff_rk2 = crk2(:,:,:,end) - csol(:,:,:);
        surf(diff_rk2(:,:,height,end)); title('RK2 Error'); az = 0; el = 90; view(az, el);
        drawnow
        colorbar;
        subplot(3,2,6)
        diff_c2l = c2tl(:,:,:,end) - csol(:,:,:);
        surf(diff_c2l(:,:,height,end)); title('2TL Error'); az = 0; el = 90; view(az, el);
        drawnow
        colorbar;
        %pause
    end
    dx = dx/2;
    dt = dt/2;
end % 
end
end

%save([RES_PATH,'error-num-rk2-',datestr(now)], 'error_rk2');

% plot the error
m = {'h','o','*','.','x','s','d','^','v','>','<','p','h'};
ls = {'-','--', ':', '-.'};
cmap = {'r','g', 'b', 'c', 'm', 'y', 'k' , 'w'};
cmapcnt = 1;
legendInfo = {};
s_fig_name = ['error_cmp_', datestr(now)];
f = figure('Name',s_fig_name,'units','normalized','outerposition',[0 0 1 1]);
for interptypecnt=1:length(INTERP_LIST)
    for fcnt = 1:length(VFREQ_LIST)
        freq = VFREQ_LIST(fcnt);
        % create and save the error graph
        figure(f)
        h = semilogy(n_level,error_rk2(interptypecnt,:,fcnt));
        mi = mod(cmapcnt,length(m)) + 1;
        set(h, 'Marker', m{mi});
        lsi = mod(cmapcnt,length(ls)) + 1;
        set(h, 'LineStyle', ls{lsi});
        ci = mod(cmapcnt, length(cmap)) + 1;
        set(h, 'Color', cmap{ci});
        cmapcnt= cmapcnt+1;
        mes = sprintf('rk2-%s-vfreq%d',INTERP_LIST{interptypecnt},freq);
        legendInfo{end+1}= mes;
        hold on;
        i = semilogy(n_level,error_2tl(interptypecnt,:,fcnt));
        mi = mod(cmapcnt,length(m)) + 1;
        set(i, 'Marker', m{mi});
        lsi = mod(cmapcnt,length(ls)) + 1;
        set(i, 'LineStyle', ls{lsi});
        ci = mod(cmapcnt, length(cmap)) + 1;
        set(i, 'Color', cmap{ci});
        cmapcnt= cmapcnt+1;
        mes = sprintf('twotl-%s-vfreq%d',INTERP_LIST{interptypecnt},freq);
        legendInfo{end+1}= mes;
        hold on;
    end
end
legend(legendInfo);
title('Error Comparision');
xlabel('Domain/Time Resolution');
ylabel('Error (logarithmic)');
s_fig_path = [RES_PATH,s_fig_name,fig_format];
saveas(f,s_fig_path);

end