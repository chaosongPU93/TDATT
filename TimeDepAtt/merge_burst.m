function [ntran,mergeind,mergetran] = merge_burst(tranhf,tranlf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ntran,mergeind,mergetran] = merge_burst(tranhf,tranlf)
% This function is to merge the tremor bursts detected by hf and lf (e.g. 2
% groups). If a burst in hf and another burst in lf have an average fraction 
% of overlap in time large than a threshold, we think the 2 bursts are
% essentially the same, and they should be merged by defining the new range
% with the min start and max end time. Otherwise, both bursts are viewed as
% different and saved into a new matrix
% 
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/09/27
% Last modified date:   2021/09/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first combine and sort them
tran = [tranhf; tranlf];
tran = sortrows(tran,[1 2 3]);

ntran = [];
k = 0;
mergeind = [];
mergetran = [];

i = 1;
while i < size(tran,1)
  j = i+1;
  if tran(i,1)~=tran(j,1)
    ntran = [ntran; tran(i,:)];
  else
%     nsameday = sum(tran(:,1) == tran(i,1))
    while j<=size(tran,1) && tran(i,1)==tran(j,1) && tran(i,3)>tran(j,2)
      j = j+1;
    end
    if j>i+1  % if there is overlap, now j is at least larger than i by 2
      % no matter how overlap is, we care only the start and end of the several bursts
      k = k+1;
      mergetran(k, :) = [tran(i,1) tran(i,2) tran(j-1,3)]; % again, j denotes the next non-overlap
      mergeind{k} = i:j-1;
      ntran = [ntran; tran(i,1) tran(i,2) tran(j-1,3)];
    else
      ntran = [ntran; tran(i,:)];
    end
  end
  i = j;
end

if i == size(tran,1)
  ntran = [ntran; tran(i,:)];
end  
  

%%%%%%%%%%%%%%%  OLD algorithm, abandoned  %%%%%%%%%%%%%%%
% k = 0;  % number of merged bursts
% mergeind = [];
% mergetran = [];
% for i = 1: size(tranhf,1)
%   
%   for j = 1: size(tranlf,1)
%     if tranhf(i,1) == tranlf(j,1) && tranhf(i,3) > tranlf(j,2)
%       overlap = tranhf(i,3) - tranlf(j,2);  % overlap in time
%     elseif tranhf(i,1) == tranlf(j,1) && tranhf(i,2) < tranlf(j,3)
%       overlap = tranlf(j,3) - tranhf(i,2);  % overlap in time
%     else
%       overlap = [];
%     end
%     lenhf = tranhf(i,3) - tranhf(i,2);  % length in time 
%     lenlf = tranlf(j,3) - tranlf(j,2);
%     if ~isempty(overlap)
%       ovlapperc = (overlap/lenhf + overlap/lenlf)/2;  % fraction of overlap
%       if ovlapperc >= ovlapthr  %then we think they are essentially the same, no need to separate
%         k = k+1;
%         mergeind(k,:) = [i j];  % note the indice of bursts have overlapping pairs
%         mergetran(k,:) = [tranhf(i,1) min(tranhf(i,2),tranlf(j,2)) ...
%           max(tranhf(i,3),tranlf(j,3))];  % merged burst
%       end
%     end
%     
%   end
% end
% ntranhf = tranhf;
% ntranlf = tranlf;
% if ~isempty(mergeind)
%     ntranhf(unique(mergeind(:,1)),:) = [];  % delete the original overlapping bursts
%     ntranlf(unique(mergeind(:,2)),:) = [];
% end
% ntran = [ntranhf; ntranlf; mergetran];  % add in the merged ones
% ntran = sortrows(ntran,[1 2]);  
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
