function [u,v,w] = vel_rot(t,x,y,z,xc,yc,zc)
    global OM;
    global GVFREQ;

    tfactor = 1+sin(2*pi*GVFREQ*t);

    u = -OM*(y-xc)*tfactor;
    v =  OM*(x-yc)*tfactor;
    w = zeros(size(z));
end
