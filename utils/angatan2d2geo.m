function anggeo = angatan2d2geo(angatan2d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anggeo = angatan2d2geo(angatan2d)
%
% A transformation from the angle under the convention of 'atan2d', where 
% atan2d(Y,X) is the four quadrant arctangent of the elements of X and Y
% such that -180 <= atan2d(Y,X) <= 180, to the angle under the convention 
% of geographical north as the start and count clockwise.
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/10/05
% Last modified date:   2022/10/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temp = -angatan2d;  % swap the sign 
if temp < 0
  temp = temp+360;  % force the range be from -180,180 to 0,360, clockwise
end
temp = temp+90; % change the start direction from east to north
if temp > 360
  temp = temp-360;  % in case of overflow
end
  
anggeo = temp;