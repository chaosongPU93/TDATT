function iran = findwhichrange(data,ranges)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iran = findwhichrange(data,ranges)
%
% Imagine you have some 'data' that belong to somewhere of the 2-column 
% matrix 'ranges' in which each row represent the start and end of some data
% range (has to be non-overlapping with each other in case of ambiguity). And
% you just want to know which range this 'data' falls in.
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/05/24
% Last modified date:   2022/05/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nran = size(ranges,1);
nd = length(data);

iran = zeros(nd,1);

if nran>1 && ranges(2,1) > ranges(1,2)    % meaning that there is a gap between each sub-range  
  for id = 1: nd
    for i = 1: nran
      if data(id) >= ranges(i,1) && data(id) <= ranges(i,2)
        iran(id) = i;
        break
      else
        continue
      end
    end
    %if no match is found for some data, it will show '-inf' or 'inf' depending on which side it is
    %on the available ranges
    if data(id) < ranges(1,1)
      iran(id) = -inf;
    elseif data(id) > ranges(end,2)
      iran(id) = inf;
    end
  end
 
else    % meaning that there is no gap in between
  for id = 1: nd
    for i = 1: nran
      if data(id) >= ranges(i,1) && data(id) < ranges(i,2)
        iran(id) = i;
        break
      else
        continue
      end
    end
    if data(id) == ranges(end,2)
      iran(id) = nran;
    end
    %if no match is found for some data, it will show '-inf' or 'inf' depending on which side it is
    %on the available ranges
    if data(id) < ranges(1,1)
      iran(id) = -inf;
    elseif data(id) > ranges(end,2)
      iran(id) = inf;
    end
  end
  
end




