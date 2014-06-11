function [ output_args ] = main( input_args )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
clear all; clear globals; % constants  and preamble
addpath('./semilag/');
addpath('./tree/');
addpath('./common/');

global verbose;
global dim; dim = 2;

% RUN PARAMETERS
maxErrorPerNode = 0.1;      % Error per box
maxLevel        = 20;       % maximum tree depth
resPerNode      = 15;       % Resolution per Node
verbose         = false;

% MAIN SCRIPT
fconc = @func1;
% fvelx = @func2;
% fvely = @func3;
% fvelz = @func4;

% PLOT FUNCTION
res = 1000;
x = linspace(0,1,res);
y = linspace(0,1,res);
[xx,yy] = meshgrid(x,y);
Z1 = fconc(xx,yy); 

contour(xx,yy,Z1);
colorbar;
axis off;
hold on;

% CONCENTRATION TREE
c = qtree;
c.insert_function(fconc,maxErrorPerNode,maxLevel,resPerNode);
c.plottree

cleaves = c.leaves();
for lvcnt = 1:length(cleaves)
    cleaf = cleaves{lvcnt};
end

    function value = func1(x,y)
        xc = 0.5;
        yc = 0.5;
        value = gaussian(x,y,xc,yc);
    end
end


