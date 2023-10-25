function newcat = refine_catalog_by_trange(oldcat,daycol,seccol,trange)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to refine the old catalog by selecting only the detections
% that are inside the time ranges 'trange' whose 1st col is date (yyyyddd),
% 2nd col and 3rd is start and end time in sec.
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/08
% Last modified date:   2021/11/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
newcat = [];
for i = 1: size(trange, 1)
  temp = oldcat(oldcat(:,daycol)==trange(i,1) & oldcat(:,seccol)>=trange(i,2) & ...
                oldcat(:,seccol)<=trange(i,3), :);
  newcat = [newcat; temp];
end