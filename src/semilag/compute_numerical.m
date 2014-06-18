function [ cnumsol ] = compute_numerical(cinit, xx, yy, zz, u, v, w, ti, dt, tn, INTERP_TYPE, NUM_SCHEME)
%COMPUTE_NUMERICAL 
% 
warning('off','all');
global dim;
V2PREVTSTEP = 1;
VPREVTSTEP  = 2;
VCURTSTEP   = 3;
VNEXTSTEP   = 4;
cnumsol = cinit;
t = [ti-2*dt, ti-dt, ti, ti+dt];
fprintf('computing numerical solution:\n');
for tstep=1:tn
    fprintf('scheme: %s interp.: %s timestep: %d time:%f\n', NUM_SCHEME, INTERP_TYPE, tstep,tstep*dt);
    switch NUM_SCHEME
        case 'rk2'
            %cnumsol(:,:,:,tstep+1) = advect_sl_rk2( cnumsol, xx, yy, zz, u, v, w, t, tstep, dt, INTERP_TYPE, INTERP_TYPE );
            cnumsol(:,:,:,tstep+1) = semilag_rk2( @interp_conc,@interp_vel, xx, yy, zz, t, tstep);
        case '2tl'
            cnumsol(:,:,:,tstep+1) = semilag_2tl( @interp_conc, xx, yy, zz, u, v, w, tstep, dt, INTERP_TYPE, INTERP_TYPE);
        otherwise
            error('Numerical scheme is unknown.');
    end
    assert(sum(sum(sum(isnan(cnumsol(:,:,:,tstep+1))))) == 0,'NaN found in solution data.');

    % clear the precomputed velociy values for the previous time step
    interp_vel(t,xx,yy,zz, true);
    
    % get the velocity values for the next time step
    t = t + dt;
    for tcnt=1:length(t)
        [u(:,:,:,tcnt), v(:,:,:,tcnt), w(:,:,:,tcnt)] = vel_rot(t(tcnt),xx,yy,zz,0.5,0.5,0.5);
    end

end

%/* ************************************************** */
    function [ci] = interp_conc(xt,yt,zt,tstep)
        cp = cnumsol(:,:,:,tstep);
        switch dim
            case 3, ci = interp3(xx,yy,zz,cp,xt,yt,zt, INTERP_TYPE);
            case 2, ci = interp2(xx,yy,cp,xt,yt, INTERP_TYPE);
        end
        % TODO: fix this. gaussian function at time t is needed.
        % ce = gaussian( xt, yt, zt, 0.5, 0.5, 0.5, 0);
        out = xt<0 | xt>1  | yt<0 | yt>1 | zt<0 | zt>1;
        ci(out) = 0;%ce(out);
    end

%/* ************************************************** */
    function [uq,vq,wq] = interp_vel(tq,xq,yq,zq,clear_values)
        if nargin < 5, clear_values = false; end;                
        persistent ut vt wt;
        if clear_values, clear ut vt wt; return; end;
        
        % precomputing the interpolated velocity values
        n = 10;
        tau = (t(VCURTSTEP)-t(VPREVTSTEP))/n;
        tqg = linspace(t(VPREVTSTEP),t(VCURTSTEP),2*n+1);        
        if (isempty(ut) || isempty(vt) || isempty(wt))
            fprintf('precomputing the interolated velocity values ... ');
            [ut, vt, wt] = interp_vel_temporal(u,v,w,t,tqg);
            fprintf('done\n');
        end
        
        iind = time2index(tq);
        [uq, vq, wq]  = interp_vel_spatial(xx,yy,zz,ut(:,:,:,iind),vt(:,:,:,iind),wt(:,:,:,iind),xq,yq,zq,INTERP_TYPE);
        
        %[ut, vt, wt] = interp_vel_temporal(u,v,w,t,tq);
        %[uq, vq, wq]  = interp_vel_spatial(xx,yy,zz,ut,vt,wt,xq,yq,zq,INTERP_TYPE);
        
        %/* ************************************************** */
        function [taui] = time2index(tq)
            ivel_indices = find(abs(tqg - tq) < 0.25*tau);
            taui = ivel_indices(1);
        end
        
        if clear_values, clear ut vt wt; end;
    end

%/* ************************************************** */
    function [uq, vq, wq] = interp_vel_spatial(xx,yy,zz,ur,vr,wr,xq,yq,zq,INTERP_TYPE)
        %INTERP_VEL_SPATIAL Spatial interpolation of velocity values.
        %
        % ur, vr, wr are the velocity values at the grid points
        % xq, yq, zq are the positions of queried velocities.
        
        switch dim
            case 3
                uq = interp3(xx,yy,zz,ur,xq,yq,zq, INTERP_TYPE);
                vq = interp3(xx,yy,zz,vr,xq,yq,zq, INTERP_TYPE);
                wq = interp3(xx,yy,zz,wr,xq,yq,zq, INTERP_TYPE);
            case 2
                uq = interp2(xx,yy,ur,xq,yq, INTERP_TYPE);
                vq = interp2(xx,yy,vr,xq,yq, INTERP_TYPE);
                wq = interp2(xx,yy,wr,xq,yq, INTERP_TYPE);
        end
        
        out = xq<0 | xq>1  | yq<0 | yq>1 | zq<0 | zq>1;
        % TODO
        [ue, ve, we] = vel_rot(0,xq,yq,zq,0.5,0.5,0.5);
        uq(out) = -ue(out);
        vq(out) = -ve(out);
        wq(out) = -we(out);
    end

%/* ************************************************** */
    function [ut, vt, wt] = interp_vel_temporal(u,v,w,t,tq)
        %INTERP_VEL_TEMPORAL Temporal interpolation of velocity values.
        %
        % tq is the time that velocity values are queried.
        % TODO: vectorize the 1D interpolation
        xsize = size(u,1);
        ysize = size(u,2);
        zsize = size(u,3);
        tsize = size(u,4);
        ut = zeros(xsize, ysize, zsize, size(tq,2));
        vt = ut;
        wt = ut;
        for i=1:xsize
            for j=1:ysize
                for k=1:zsize
                    ulist = zeros(size(t));
                    vlist = ulist;
                    wlist = ulist;
                    for cnt=1:tsize
                        ulist(cnt) =u(i,j,k,cnt);
                        vlist(cnt) =v(i,j,k,cnt);
                        wlist(cnt) =w(i,j,k,cnt);
                    end
                    ut(i,j,k,:) = interp1(t,ulist,tq,INTERP_TYPE);
                    vt(i,j,k,:) = interp1(t,vlist,tq,INTERP_TYPE);
                    wt(i,j,k,:) = interp1(t,wlist,tq,INTERP_TYPE);
                end
            end
        end
    end
end

