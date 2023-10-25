function [x, y] = ellipse_chao(xc,yc,semia,semib,dx,rotang,rotcnt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [x, y] = ellipse_chao(x0,y0,semia,semib,dx,rotang,rotcnt)
%
% This function is to obtain the points on the perimeter of a ellipse with a
% semi-major axis of 'semia' and semi-minor axis of 'semib', centered at 
% (xc,yc). The spacing of x coordinate is set as 'dx'. Also has the option
% to rotate counter-clockwise by 'rotang' degrees relative to 
% a center 'rotcnt'.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/14
% Last modified date:   2022/03/14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('rotang',0); %Set the default to no rotation
defval('rotcnt',[xc,yc]); %Set the default to rotation wrt the center

%the positive half
xup = -semia: dx: semia;
yup = semib*sqrt(1-(xup/semia).^2);

%include the other half
x = reshape([xup,fliplr(xup(1:end-1))], [], 1);
y = reshape([yup,-fliplr(yup(1:end-1))], [], 1);

%shift the center from 0,0 to real center xc,yc
x = x+xc;
y = y+yc;

%if needs to rotate counter-clockwiseby 'angrot' degrees, relative to 'rotcnt' 
[x,y] = complex_rot(x,y,rotang,rotcnt);


