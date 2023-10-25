function avgtrace = Runningavg(trace,win)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avgtrace = Runningavg(trace,win)
% This is a simple function to carry out the running average of some trace,
% single vector and matrix using a specified window. If the trace is a matrix,
% then the operation is done to each column. The window is a weight vector 
% of any type, rectangle if you want just an arithmetric mean, Gaussian, etc.
% But it should not be the same length as the trace. The kernel operation
% is a convolution between the trace and win, so that each point of 'avgtrace' 
% is a dot product between 'trace' and 'win' (weight).
% 
% 
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/03/24
% Last modified date:   2022/03/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avgtrace = zeros(size(trace));
for i = 1: size(trace,2)
  avgtrace(:,i) = conv(trace(:,i),win,'same');
end