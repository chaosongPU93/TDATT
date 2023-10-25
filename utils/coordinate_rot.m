function [newx,newy] = coordinate_rot(oldx,oldy,angle,rotcnt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [newx,newy] = coordinate_rot(oldx,oldy,angle,rotcnt)
%
% This function is used to rotate the 2-D coodinate system X-Y counter-clockwise 
% by an 'angle' to get the location of points (oldx,oldy) in the new coord system 
% X'-Y'. A shift of the center (sftx,sfty), ie., a translation is optional.
% For 3-D coodinate rotation, see 'coordinate_rot3d.m', or considering rotating
% the 2-D system multiple times about the other axis which remains unchanged every
% time.
% Note the difference between 'complex_rot'.  
%
% INPUT:  
%   oldx:     original x, one point or vector
%   oldy:     original y, one point or vector
%   angle:  calculated counter-clockwise. in degree, add '-' to incate clockwise rotation
%   rotcnt:   rotation center in the original X-Y-Z, default is (0,0)
% 
% OUTPUT:
%   newx:   rotated x in the new coordinate system
%   newy:   rotated y in the new coordinate system
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/10/02
% Last modified date:   2022/03/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%some defaults
defval('rotcnt',[0,0]);

%change degree to radian
angrad = deg2rad(angle);

%rotation matrix
R = [cos(angrad) sin(angrad); -sin(angrad) cos(angrad)];   

oldx = reshape(oldx-rotcnt(1),1,[]);  % shift to the rotation center, and force to colomun vector
oldy = reshape(oldy-rotcnt(2),1,[]);

newcoord = R * [oldx; oldy];    % matrix multiplication

newx = newcoord(1,:)';
newy = newcoord(2,:)';


