function [ c ] = slotted_cylinder( xi, xf, xx, yy, zz, dx)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
domain_length = xf - xi;
xcenter = (xf + xi) / 2;
ycenter = xcenter;
x = xi:dx:xf;
s = length(x);
radius = 0.3*domain_length;
c = zeros(size(xx));
% cylinder
xd = xx - xcenter;
yd = yy - ycenter;
dist = sqrt(xd.*xd + yd.*yd);
cmask = dist < radius;
c(cmask) = 1;
% slotted part
slottDX = radius/2;
sxi = xi + 0.5*domain_length - 0.5*slottDX;
sxf = xi + 0.5*domain_length + 0.5*slottDX;
syi = xi;
syf = xi + 2*domain_length/3;
smask = xx > sxi & xx < sxf & yy > syi & yy < syf;
c(smask) = 0;
end

