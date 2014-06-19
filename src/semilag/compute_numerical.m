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
            cnumsol(:,:,:,tstep+1) = semilag_rk2(@interp_conc,@interp_vel, xx, yy, zz, t,tstep);
        case '2tl'
            cnumsol(:,:,:,tstep+1) = semilag_2tl(@interp_conc,xx,yy,zz,u,v,w,tstep,dt,INTERP_TYPE,INTERP_TYPE);
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
    function [ci] = interp_conc(tstep,xt,yt,zt)        
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
            [ut, vt, wt] = interp_vel_temporal(u,v,w,t,tqg,INTERP_TYPE);
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
    end
end

