function [newx,newy,newz] = coordinate_rot3d(oldx,oldy,oldz,angle,rotaxis,rotcnt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to rotate the 3-D coodinate system X-Y-Z about one of 
% the 3 axes counter-clockwise by an 'angle' to get the location of points 
% (oldx,oldy) in the new coord system X'-Y'-Z'. A shift of the center 
% 'cntsft', ie., a translation is optional.
% For 2-D coodinate rotation, see 'coordinate_rot.m' 
%
% INPUT:  
%   oldx:     original x, one point or vector
%   oldy:     original y, one point or vector
%   oldz:     original z, one point or vector
%   angle:    calculated counter-clockwise. in degree, add '-' to incate clockwise rotation
%   rotaxis:  which axis the rotation is about? 1 for X, 2 for Y, 3 for Z (default) 
%   rotcnt:   rotation center in the original X-Y-Z, default is (0,0,0)
% 
% OUTPUT:
%   newx:   rotated x in the new coordinate system
%   newy:   rotated y in the new coordinate system
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/17
% Last modified date:   2022/03/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%some defaults
defval('rotaxis',3);  % default rotation is about Z axis
defval('rotcnt',zeros(3,1));  % default rotation center is (0,0,0)

%change degree to radian
angrad = deg2rad(angle);

%determine rotation matrix based on which axis the rotation is about, meaning that coordinate does
%not change
switch rotaxis
  case 1
    R = [1 0 0; 0 cos(angrad) sin(angrad); 0 -sin(angrad) cos(angrad)];
  case 2
    R = [cos(angrad) 0 sin(angrad); 0 1 0; -sin(angrad) 0 cos(angrad)];
  case 3
    R = [cos(angrad) sin(angrad) 0; -sin(angrad) cos(angrad) 0; 0 0 1];  
end

oldx = reshape(oldx-rotcnt(1),1,[]);  % shift to the rotation center, and force to colomun vector
oldy = reshape(oldy-rotcnt(2),1,[]);
oldz = reshape(oldz-rotcnt(3),1,[]);

newcoord = R * [oldx; oldy; oldz];    % matrix multiplication

newx = newcoord(1,:)';
newy = newcoord(2,:)';
newz = newcoord(3,:)';

