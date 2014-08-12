function [ u, v, w, c ] = init_fields(xi, xf, dx, xx, yy, zz, t, VF_TYPE, CF_TYPE)
xc = (xf + xi) / 2;
yc = xc;
zc = xc;
c = zeros(size(xx));

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
            case 0, c = slotted_cylinder( xi, xf, xx, yy, zz);
        end
    end
end

