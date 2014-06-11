function [ cnew ] = advect_sl_rk2( c, xx, yy, zz, u, v, w, t, tstep,dt, vel_interp_method, conc_interp_method )
%ADVECT_SL_RK2 Advect the c values one time step by using semi-lagrangian scheme
%
% The trajectory in semi-lagrangian is computed with second order runge
% kutta scheme
warning('off','all');
global dim;

%warning;
%conc_interp_method = 'spline';
V2PREVTSTEP = 1;
VPREVTSTEP  = 2;
VCURTSTEP   = 3;
VNEXTSTEP   = 4;
xsize = size(u,1);
ysize = size(u,2);
zsize = size(u,3);
tsize = size(t,2);
uback = -u;
vback = -v;
wback = -w;
n = 10;
taui = 0;
tauf = dt;
tau = (tauf - taui)/n; 
xtau = xx;
ytau = yy;
ztau = zz;
utau = uback(:,:,:,VCURTSTEP);
vtau = vback(:,:,:,VCURTSTEP);
wtau = wback(:,:,:,VCURTSTEP);
tcur = t(VCURTSTEP);

% precomputing the interpolated velocity values
fprintf('precomputing the interolated velocity values ... ');
tqg = linspace(t(VPREVTSTEP),t(VCURTSTEP),2*n+1);
[uqg, vqg, wqg] = interp_vel_temporal(tqg);
fprintf('done\n');

%figure
for taustep=1:n 
    [xtau, ytau, ztau] = rk2(xtau, ytau, ztau, utau, vtau, wtau, tau);
    tcur = tcur - tau;
    [utau, vtau, wtau] = interp_vel(xtau, ytau, ztau, tcur);
%     MS='MarkerSize';
%     plot(xx(:),yy(:),'o',MS,5); hold on;
%     plot(xtau(:),ytau(:),'gs',MS,8);  axis off; axis equal;
%     title('Traj. Points of RK2');
%     hold on
%     pause
end

% hold on
% MS='MarkerSize';
% plot(xx(:),yy(:),'o',MS,5); hold on;
% plot(xtau(:),ytau(:),'gs',MS,8);  axis off; axis equal;
% title('Traj. Points of RK2');
% hold on
% pause

cnew = interp_conc(c(:,:,:,tstep),xtau,ytau,ztau);

    function [taui] = time2index(tq)
        ivel_indices = find(abs(tqg - tq) < 0.25*tau);
        taui = ivel_indices(1);
    end

    function [xt,yt,zt] = rk2(x,y,z,u1,v1,w1,ddt)
        % dx/dt = vel(x,t)
        % k1    = dt * vel(x_n, t)
        % k2    = dt * vel(x_n+k1/2, t_n+dt/2)
        % x_n+1 = x_n + k2
        k1x = ddt*u1;
        k1y = ddt*v1;
        k1z = ddt*w1;
        % interpolate velocity in time: t_n+dt/2
        % interpolate velocity in space: x_n+k1/2
        [uts2,vts2,wts2] = interp_vel(x+0.5*k1x, y+0.5*k1y, z+0.5*k1z, tcur - 0.5*ddt);
        xt = x + ddt*uts2;
        yt = y + ddt*vts2;
        zt = z + ddt*wts2;
    end

    function [ui,vi,wi] = interp_vel(xq,yq,zq,tq)
        %[uti, vti, zti] = interp_vel_temporal(tq);
        iind = time2index(tq);
        [ui, vi, wi]  = interp_vel_spatial(uqg(:,:,:,iind), vqg(:,:,:,iind), wqg(:,:,:,iind), xq, yq, zq);
    end

    function [us, vs, ws] = interp_vel_spatial(ur, vr, wr, xq, yq, zq)
        %INTERP_VEL_SPATIAL Spatial interpolation of velocity values. 
        % 
        % ur, vr, wr are the velocity values at the grid points
        % xq, yq, zq are the positions of queried velocities.
        switch dim
            case 3
                us = interp3(xx,yy,zz,ur,xq,yq,zq, vel_interp_method);
                vs = interp3(xx,yy,zz,vr,xq,yq,zq, vel_interp_method);
                ws = interp3(xx,yy,zz,wr,xq,yq,zq, vel_interp_method);
            case 2
                us = interp2(xx,yy,ur,xq,yq, vel_interp_method);
                vs = interp2(xx,yy,vr,xq,yq, vel_interp_method);
                ws = interp2(xx,yy,wr,xq,yq, vel_interp_method);
        end
        out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
        [ue, ve, we] = vel_rot(tcur,xq,yq,zq,0.5,0.5,0.5);
        us(out) = -ue(out);
        vs(out) = -ve(out);
        ws(out) = -we(out);
    end

    function [ut, vt, wt] = interp_vel_temporal(tq)
        %INTERP_VEL_TEMPORAL Temporal interpolation of velocity values.
        %
        % tq is the time that velocity values are queried.
        % TODO: vectorize the 1D interpolation
        ut = zeros(xsize, ysize, zsize, size(tq,2));
        vt = ut;
        wt = ut;
        for i=1:xsize
            for j=1:ysize
                for k=1:zsize
                    ulist = zeros(size(t));
                    vlist = ulist;
                    wlist = ulist;
                    for tcnt=1:tsize
                        ulist(tcnt) =uback(i,j,k,tcnt);
                        vlist(tcnt) =vback(i,j,k,tcnt);
                        wlist(tcnt) =wback(i,j,k,tcnt);
                    end
                    ut(i,j,k,:) = interp1(t,ulist,tq,'cubic');
                    vt(i,j,k,:) = interp1(t,vlist,tq,'cubic');
                    wt(i,j,k,:) = interp1(t,wlist,tq,'cubic');
                end
            end
        end
    end

    function [ci] = interp_conc(cp,xt,yt,zt)
        switch dim
            case 3, ci = interp3(xx,yy,zz,cp,xt,yt,zt, conc_interp_method);
            case 2, ci = interp2(xx,yy,cp,xt,yt, conc_interp_method);
        end
        % TODO: fix this. gaussian function at time t is needed.
        % ce = gaussian( xt, yt, zt, 0.5, 0.5, 0.5, 0);
        out = xt<0 | xt>1  | yt<0 | yt>1 | zt<0 | zt>1;
        ci(out) = 0;%ce(out);
    end
end
