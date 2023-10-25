function [xascend,ikeep,idiscard] = ascendvector(x,wgt,wgtthr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [xchron,ikeep,idiscard] = ascendvector(x,wgt)
%
% Imgine there is a data set has a general ascending order, but there are 
% some points deviating from this general trend. Now you want to discard 
% these points to enforce the ascending order. The default is each data
% point has an equal weight of 1, but you can also assign a weight vector
% so that points which deviate from the ascending order but has a weight
% higher than some threshold would still be preserved.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/02/06
% Last modified date:   2023/02/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('wgt',ones(size(x,1)));
defval('wgtthr',0);

npts = size(x,1);

idiscard = [];
iref = 1;
ikeep = [iref];

if ~isequaln(wgt,ones(size(x,1)))
  %if all data points have unequal weights
  while iref <= npts-1 
    ichk = iref+1;
    if x(ichk) >= x(iref)
      ikeep = [ikeep; ichk];
      iref = ichk;
    else
      while x(ichk) < x(iref)
        if wgt(ichk) < wgtthr
          idiscard = [idiscard; ichk];
          ichk = ichk+1;
          if ichk > npts
            break
          end
        else
          ikeep = [ikeep; ichk];
          iref = ichk;
          break
        end
      end
      ikeep = [ikeep; ichk];
      iref = ichk;
    end
  end
  
else
  %if all data points have equal weight of 1
  while iref <= npts-1
    ichk = iref+1;
    if x(ichk) >= x(iref)
      ikeep = [ikeep; ichk];
      iref = ichk;
    else
      while x(ichk) < x(iref)
        idiscard = [idiscard; ichk];
        ichk = ichk+1;
        if ichk > npts
          break
        end
      end
      ikeep = [ikeep; ichk];
      iref = ichk;
    end
  end
  
end
xascend = x(ikeep);



