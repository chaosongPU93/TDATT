function M0=mag2moment(Mw)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M0=mag2moment(Mw)
% This script converts moment magnitude to scalar moment in unit 'Nm'
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/08/08
% Last modified date:   2023/08/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%using equation below
M0 = 10.^(1.5*(10.7+Mw))/1.0e7;
