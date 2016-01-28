function [u,v,w] = vel_rot(t,x,y,z,xc,yc,zc)
global om;
global gvfreq;

tfactor = 1+sin(2*pi*gvfreq*t);

u = -om*(y-xc)*tfactor;
v =  om*(x-yc)*tfactor;
w = zeros(size(z));
end
