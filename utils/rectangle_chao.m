function [x, y] = rectangle_chao(xc,yc,wid,hgt,inc,rotang,rotcnt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x, y] = rectangle_chao(x0,y0,semia,semib,dx)
%
% This function is to obtain the points on the perimeter of a rectangle with a
% width of 'wid' and height of 'hgt', centered at (xc,yc).
% The spacing of each coordinate is set as 'inc'. Also has the option
% to rotate counter-clockwise by 'rotang' degrees relative to a center 'rotcnt'.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/14
% Last modified date:   2022/03/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('rotang',0); % set the default to no rotation
defval('rotcnt',[0 0]); % default rotation center is (0,0)

%Top 
x1 = (-wid/2: inc: wid/2)';
y1 = hgt/2*ones(length(x1),1);

%Right
y2 = (hgt/2-inc: -inc: -hgt/2+inc)';
x2 = wid/2*ones(length(y2),1);

%Bottom
x3 = flipud(x1);
y3 = -y1;

%Left
y4 = flipud(y2);
x4 = -x2;

%Assembly
x = [x1; x2; x3; x4; x1(1)];
y = [y1; y2; y3; y4; y1(1)];

%shift the center from 0,0 to real center xc,yc
x = x+xc;
y = y+yc;

%if needs to rotate counter-clockwiseby 'angrot' degrees, relative to 'rotcnt' 
[x,y] = complex_rot(x,y,rotang,rotcnt);




