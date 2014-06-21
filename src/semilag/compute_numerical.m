function [ cnumsol ] = compute_numerical(cinit, xx, yy, zz, u, v, w, ti, dt, tn, INTERP_TYPE, NUM_SCHEME)
%COMPUTE_NUMERICAL 
% 
warning('off','all');
% V2PREVTSTEP = 1;
% VPREVTSTEP  = 2;
% VCURTSTEP   = 3;
% VNEXTSTEP   = 4;
t = [ti-2*dt, ti-dt, ti, ti+dt];
cnumsol = cinit;

fprintf('computing numerical solution:\n');
for tstep=1:tn
    fprintf('scheme: %s interp.: %s timestep: %d time:%f\n', NUM_SCHEME, INTERP_TYPE, tstep,tstep*dt);
    switch NUM_SCHEME
        case 'rk2'
            cnumsol(:,:,:,tstep+1) = semilag_rk2(cnumsol,xx,yy,zz,@interp_conc,u,v,w,t,@interp_vel_precomputed,tstep,INTERP_TYPE);
            % clear the precomputed velociy values for the previous time step
            interp_vel_precomputed(xx,yy,zz,u,v,w,t,0,0,0,0,INTERP_TYPE,true);
        case '2tl'
            cnumsol(:,:,:,tstep+1) = semilag_2tl(cnumsol,xx,yy,zz,@interp_conc,u,v,w,tstep,dt,INTERP_TYPE);
        otherwise
            error('Numerical scheme is unknown.');
    end
    assert(sum(sum(sum(isnan(cnumsol(:,:,:,tstep+1))))) == 0,'NaN found in solution data.');
    % get the velocity values for the next time step
    t = t + dt;
    for tcnt=1:length(t)
        [u(:,:,:,tcnt), v(:,:,:,tcnt), w(:,:,:,tcnt)] = vel_rot(t(tcnt),xx,yy,zz,0.5,0.5,0.5);
    end
end
end

