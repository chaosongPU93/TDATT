function [xrot,yrot] = complex_rot(x,y,angle,rotcnt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xrot,yrot] = complex_rot(x,y,angle,rotcnt)
%
% This function is used to rotate the original point by an angle counter-
% clockwise relative to the rotation center, and output the coordinates 
% in the original system. No matter how the rotation is, they
% are all in the same coordinate system.
% Note the difference between 'coordinate_rot'.  
%
% INPUT:  
%   x:     original x
%   y:     original y
%   angle:  calculated counter-clockwise. in degree, add '-' to incate clockwise rotation
%   rotcnt:    rotation center (xrc,yrc)
%   
% 
% OUTPUT:
%   xrot:   rotated x in the same coordinate system
%   yrot:   rotated y in the same coordinate system
%
% Chao Song, chaosong@princeton.edu
% First created date:   2019/10/02
% Last modified date:   2019/10/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('rotcnt',[0 0]); % default rotation center is (0,0)

% c = x0+1i.*y0;
x = x-rotcnt(1);  % x relative to rotation center
y = y-rotcnt(2);
c = complex(x,y);
crot = c.*exp(1i * deg2rad(angle));   % rotate counter-clockwise by angle
xrot = real(crot);
yrot = imag(crot);

%revert to original coordinate from that relative to rotation center
xrot = xrot+rotcnt(1);  % x relative to rotation center
yrot = yrot+rotcnt(2);
