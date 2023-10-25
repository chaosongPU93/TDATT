function stats = statsofrobustlnfit(fitobj,gof,output,x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stats = statsofrobustlnfit(fitobj,gof,output,x,y)
%
% This is to obtain some more statistics of the robust linear fitting of data
% set x and y and returned objects.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/26
% Last modified date:   2023/01/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Some statistics
%compute the HF weights in robust linear regression
% x = tsplst/sps;
res = output.residuals;   % usual residuals
hatmat = x/inv(x'*x)*x';
h = zeros(size(hatmat,1),1);    % leverage of least square
for jj = 1 : size(hatmat,1)
  h(jj) = hatmat(jj,jj);
end
radj = res./sqrt(1-h);      % adjusted residuals
K = 4.685;
s = mad(res,1)/0.6745;
u = radj/(K*s);
wt = zeros(length(u),1);    % rubust weight of next iteration
for jj = 1 : length(u)
  if abs(u(jj)) < 1
    wt(jj) = (1-(u(jj))^2)^2;
  else
    wt(jj) = 0;
  end
end
%get the standard error of the estimated parameters, may indicate the compare the quality
%of fitting cross different RTMs, NOTE that rmse, mse are definitely not a good measure
%the obtained CI is comfirmed to be correct by comparing it with the 'fitobjhfprop'
coef = coeffvalues(fitobj);
slope = coef(1);
intcpt = coef(2);
slopese = gof.rmse./sqrt(sum((x-mean(x)).^2));
est = slope;
slopeCI = confidence_interval_general(est,slopese,length(x)-2,95);
intcptse = slopese.*sqrt(sum(x.^2)./length(x));
est = intcpt;
intcptCI = confidence_interval_general(est,intcptse,length(x)-2,95);

%compute weighted pearson coeff, which is to evaluate the linear correlation between x and y,
%algorithm here is related to built-in "corr(x,y,'type','Pearson')"
x_bar = wt_mean(x,wt);
y_bar = wt_mean(y,wt);
x_var = sum(wt.*(x-x_bar).^2) / sum(wt);
y_var = sum(wt.*(y-y_bar).^2) / sum(wt);
xy_cov = sum(wt.*(x-x_bar).*(y-y_bar)) / sum(wt);
pearwt = xy_cov / sqrt(x_var*y_var);

%save statistics into a structure
stats.slope = slope;
stats.slopese = slopese;
stats.slopeCI = slopeCI;
stats.intcpt = intcpt;
stats.intcptse = intcptse;
stats.intcptpropCI = intcptCI;
stats.wt = wt;
stats.pearwt = pearwt;
stats.fitobj = fitobj;
stats.gof = gof;
stats.output = output;


