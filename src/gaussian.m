function [ conc ] = gaussian(X, Y, xcenter, ycenter, theta, sigma_x, sigma_y)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<3, xcenter = 0.5; end;
if nargin<4, ycenter = 0.5; end;
if nargin<5, theta=0; end;
if nargin<6, sigma_x = 0.05; end;
if nargin<7, sigma_y = 0.05; end;

A = 1;
x0 = xcenter; 
y0 = ycenter;
a = cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
c = sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;
conc = A*exp( - (a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) ;
end

