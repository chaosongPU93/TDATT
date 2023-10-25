function sumvect = sumofrange(vect,nsum,direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sumvect = sumofrange(vect,nsum,direction)
%
% Sometimes when you have a vector 'vect', and you want to know the sum of 
% 'nsum' consecutive elements along an custom 'direction', 'sumofrange' aims
% to solve such a problem.
% 'direction' can be 'center' (default; equal number forward and backward), 
% or 'forward' (count nsum from the current element forward), or 'backward'
% (count nsum from the current element backward).
% The output 'sumvect' has the same length as the input 'vect', but note some
% elements in 'sumvect' does not have the same number of 'nsum' contributing
% to the sum, depending on the 'direction'.
%
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/09
% Last modified date:   2022/06/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('direction','center');

N = length(vect);
sumvect = zeros(N, 1);

if isequal(direction,'forward')
  for i = 1: N
    sumvect(i) = sum(vect(max(i-nsum+1,1): i));
  end
elseif isequal(direction,'backward')
  for i = 1: N
    sumvect(i) = sum(vect(i: min(i+nsum-1,N)));
  end
elseif isequal(direction,'center')
  if rem(nsum,2)==1
    for i = 1: N
      sumvect(i) = sum(vect(max(i-round((nsum-1)/2),1): min(i+round((nsum-1)/2),N)));
    end
  else
    for i = 1: N
      sumvect(i) = sum(vect(max(i-round(nsum/2),1): min(i+round((nsum-2)/2),N)));
    end
  end
end








