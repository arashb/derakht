function [uq, vq, wq] = interp_vel_temporal(u,v,w,t,tq,INTERP_TYPE)
%INTERP_VEL_TEMPORAL Temporal interpolation of velocity values.
%
% tq is the time that velocity values are queried.
warning('off', 'MATLAB:interp1:UsePCHIP') ;

tn = size(t,2);
tqn = size(tq,2);
s = size(u);
utmp = zeros(s(1)*s(2)*s(3),tn);
vtmp = utmp;
for cnt=1:tn
    tmp = u(:,:,:,cnt);
    utmp(:,cnt) = tmp(:);
    tmp = v(:,:,:,cnt);
    vtmp(:,cnt) = tmp(:);
end
uttmp = utmp';
vttmp = vtmp';
uqtmp = interp1(t,uttmp,tq,INTERP_TYPE);
vqtmp = interp1(t,vttmp,tq,INTERP_TYPE);

uq = zeros(s(1),s(2),s(3),tqn);
tmp = zeros(s(1),s(2),s(3));
vq = uq;
wq = uq;
for cnt=1:tqn
    tmp(:) = uqtmp(cnt,:);
    uq(:,:,:,cnt) = tmp;
    tmp(:) = vqtmp(cnt,:);
    vq(:,:,:,cnt) = tmp;
end
end
