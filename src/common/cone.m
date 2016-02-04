function [ c ] = cone( xx, yy, xc, yc, r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if nargin<3, xc = 0.65; end;
    if nargin<4, yc = 0.65; end;
    if nargin<5, r = 0.20; end;

    xd = xx - xc;
    yd = yy - yc;
    dd = xd.*xd + yd.*yd;
    c = 1 - dd./r^2;
    dist = sqrt(dd);
    cmask = dist > r;
    c(cmask) = 0;
end
