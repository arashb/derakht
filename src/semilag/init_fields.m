function [ u, v, w, c ] = init_fields(xi, xf, dx, xx, yy, zz, ti, dt, VF_TYPE, CF_TYPE)
xc = (xf + xi) / 2;
yc = xc;
zc = xc;
c = zeros(size(xx));
t = [ti-2*dt, ti-dt, ti, ti+dt];
V2PREVTSTEP = 1;
VPREVTSTEP  = 2;        % index for the vel. values of the previous time step
VCURTSTEP   = 3;        % index for the vel. values of the current time step
VNEXTSTEP   = 4;        % index for the vel. values of the next time step

% init the velocity values
init_vel();
% init the concentration values
init_conc();
assert(sum(sum(sum(isnan(c)))) == 0,'NaN found in initial concectration data.');
    
    function init_vel()
        switch VF_TYPE
            case 1
                % get the values for the needed times
                for tcnt=1:length(t)
                    [u(:,:,:,tcnt), v(:,:,:,tcnt), w(:,:,:,tcnt)] = vel_rot(t(tcnt),xx,yy,zz,xc,yc,zc);
                end
            case 0
                for tcnt=1:length(t)
                    u(:,:,:,tcnt) = ones(size(xx));
                    v(:,:,:,tcnt) = ones(size(yy));
                    w(:,:,:,tcnt) = ones(size(zz));
                end
        end
    end

    function init_conc()
        switch CF_TYPE
            case 1, c  = gaussian( xx, yy, xc, yc, 0);
            case 0, c = slotted_cylinder( xi, xf, xx, yy, zz, dx);
        end
    end
end

