function [mcoef,mlag] = xcorrmax(trace1,tracei,MAXLAG,SCALEOPT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return the value and location of max CC between trace i and trace 1 using
% 'xcorr' allowing a max lag of 'MAXLAG' with the scaling option 'SCALEOPT'. 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2023/01/03
% Last modified date:   2023/01/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defval('SCALEOPT','coeff');

[coef,lag] = xcorr(trace1,tracei,MAXLAG,SCALEOPT);
[mcoef, idx] = max(coef);   % max of master raw cc
mlag = lag(idx);   % offset in samples
