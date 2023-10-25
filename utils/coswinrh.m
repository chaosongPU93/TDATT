function w = coswinrh(len,ratio)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w = coswinrh(wlen,ratio)
%
% This function is to create an RIGHT-hand cosine tapering window. The full
% window has a length of 'len' points, and a cosine fraction of 'ratio' that 
% is tapered on the right.  See also 'coswinlh', 'tukeywin'. 
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/05/24
% Last modified date:   2022/05/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ratio >= 1 
  ratio = 1;
end

if ratio <= 0
    w = ones(len,1);
else
    t = linspace(0,1,len)';
    % Defines period of the taper as 1/2 period of a sine wave.
    th = floor(ratio*(len-1))+1;
    % Window is defined in three sections: taper, constant, taper
    w = [ ones(len-th,1); (1+cos(pi/ratio*(t(len-th+1:len) - 1 + ratio)))/2 ]; 
end