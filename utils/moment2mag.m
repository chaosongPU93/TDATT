function Mw=moment2mag(M0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mw=moment2mag(M0)
% This script converts scalar moment in unit 'Nm' to moment magnitude 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/08/08
% Last modified date:   2023/08/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mw = (2/3)*log10(M0*1.0e7)-10.7;