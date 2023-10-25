function [fitobj,gof,outp,yfit] = linear_bisquare_fit_free(xobs,yobs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to create a linear fitting object with free constraints,
% ie., invert for both slope and intercept, using the 'bisquare' robust
% algorithm, return the fit parameters, quality, and fitted y. 
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/08
% Last modified date:   2021/11/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create fit object with free constraints
fttpfree = fittype( @(a,b,x) a*x+b);

% linear robust least square
[fitobj,gof,outp] = fit(xobs, yobs,fttpfree,'Robust','Bisquare','StartPoint',[1 1]);
yfit = feval(fitobj,xobs);

