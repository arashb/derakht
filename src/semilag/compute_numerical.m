function [ cnumsol ] = compute_numerical(cinit, xx, yy, zz, u, v, w, ti, dt, tn, INTERP_TYPE, NUM_SCHEME)
%COMPUTE_NUMERICAL 
% 
cnumsol = cinit;
t = [ti-2*dt, ti-dt, ti, ti+dt];
%f = figure;
fprintf('computing numerical solution:\n');
for tstep=1:tn
    fprintf('scheme: %s interp.: %s timestep: %d time:%f\n', NUM_SCHEME, INTERP_TYPE, tstep,tstep*dt);
    switch NUM_SCHEME
        case 'rk2'
            cnumsol(:,:,:,tstep+1) = advect_sl_rk2( cnumsol, xx, yy, zz, u, v, w, t, tstep, dt, INTERP_TYPE, INTERP_TYPE );
        case '2tl'
            cnumsol(:,:,:,tstep+1) = advect_sl_2tl( cnumsol, xx, yy, zz, u, v, w, tstep, dt, INTERP_TYPE, INTERP_TYPE);
        otherwise
            error('Numerical scheme is unknown.');
    end
    assert(sum(sum(sum(isnan(cnumsol(:,:,:,tstep+1))))) == 0,'NaN found in solution data.');

    % get the velocity values for the next time step
    t = t + dt;
    % TODO: just shift the previous velocities and compute only velocity
    % for ti+dt
    for tcnt=1:length(t)
        [u(:,:,:,tcnt), v(:,:,:,tcnt), w(:,:,:,tcnt)] = vel_rot(t(tcnt),xx,yy,zz,0.5,0.5,0.5);
%         figure(f)
%         quiver(u(:,:,2,tcnt),v(:,:,2,tcnt));
%         grid on
%         pause
    end

end

end

