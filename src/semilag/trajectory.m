function [ xxorigin,yyorigin,zzorigin ] = trajectory(xx, yy, zz, velfunc, ti, tf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
xc = 0.5;
yc = xc;
zc = 0;
nx = size(xx,1);
ny = size(yy,2);
nz = size(zz,3);
p2c=@(p)reshape(p,nx,ny,nz,3);
c2p=@(x,y,z)[x(:);y(:);z(:)];
pos0 = c2p(xx,yy,zz);
fprintf('computing trajectory\n');
opt=odeset('RelTol',0.5,'AbsTol',0.5);
[T,Y] = ode45(@rhs, [tf ti], pos0,opt);
pos=p2c(Y(end,:));
xxorigin=pos(:,:,:,1);
yyorigin=pos(:,:,:,2);
zzorigin=pos(:,:,:,3);

    function r=rhs(t,pos)
        pos = p2c(pos);
        x=pos(:,:,:,1); y = pos(:,:,:,2); z = pos(:,:,:,3);
        [ru,rv,rw]= velfunc(t,x,y,z,xc,yc,zc);
        r = c2p(ru,rv,rw);
    end
end