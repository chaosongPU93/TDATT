function dfdx = dfxdx(fx,x)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is used to compute finite difference sense of derivative of
% f(x) with respect to x, at each sample. Can be also related to the 1-D 
% gradient of f(x) with respect to x, if the direction is considered.
% 
% 
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2022/03/17
% Last modified date:   2022/03/17
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fx = reshape(fx,[],1);
x = reshape(x,[],1);

dfx = diff(fx);
dx = diff(x);
dfdx = dfx./dx;
  
