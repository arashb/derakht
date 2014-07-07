function [uq,vq,wq] = interp_vel_precomputed(xx,yy,zz,u,v,w,t,tq,xq,yq,zq,INTERP_TYPE,fout,clear_values)
if nargin < 14, clear_values = false; end;
persistent ut vt wt;
if clear_values, clear ut vt wt; return; end;

V2PREVTSTEP = 1;
VPREVTSTEP  = 2;
VCURTSTEP   = 3;
VNEXTSTEP   = 4;

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
[uq, vq, wq]  = interp_vel_spatial(xx,yy,zz,ut(:,:,:,iind),vt(:,:,:,iind),wt(:,:,:,iind),xq,yq,zq,INTERP_TYPE,fout);

    %/* ************************************************** */
    function [taui] = time2index(tq)
        ivel_indices = find(abs(tqg - tq) < 0.25*tau);
        taui = ivel_indices(1);
    end
end