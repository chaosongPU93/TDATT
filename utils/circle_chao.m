function [x, y] = circle_chao(x0,y0,radius,ddeg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x, y] = circle_chao(x0,y0,radius,ddeg)
%
% This function is to obtain the points on the perimeter of a circle with a
% radius centered at (x0,y0). The spacing of angle in degree is set as drad
%
% Chao Song, chaosong@princeton.edu
% First created date:   2020/08/13
% Last modified date:   2020/08/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drad = deg2rad(ddeg);
theta = 0: drad: 2*pi;
x = radius * cos(theta) + x0;
y = radius * sin(theta) + y0;


