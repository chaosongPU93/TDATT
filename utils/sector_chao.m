function [x, y] = sector_chao(x0,y0,radius,stdeg,eddeg,ddeg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to obtain the points on the perimeter of a sector (part 
% of a circle) with a radius centered at (x0,y0) from a start angle to end
% angle. The spacing of angle in degree is set as drad. Angle is measured
% from E counterclockwise.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2020/08/19
% Last modified date:   2020/08/19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stdeg = rem(stdeg,360);
% eddeg = rem(eddeg,360);
drad = deg2rad(ddeg);
strad = deg2rad(stdeg);
edrad = deg2rad(eddeg);
theta = strad: drad: edrad;
x = radius * cos(theta) + x0;
y = radius * sin(theta) + y0;


