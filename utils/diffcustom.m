function dvect = diffcustom(vect,nsep,direction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dvect = diffcustom(vect,nsep)
%
% Sometimes you want to know the difference between the sorted elements in a
% vector 'vect', like N - (N-1), then you can use 'diff'. 'diffcustom' aims 
% to solve similar problem but for non-consecutive elements, like N-(N-2); 
% N-(N-3); N-(N-nsep)...
% You could do this along different 'direction', 'forward' [N-(N-nsep)], 
% or 'backward' [(N+nsep)-N].
% The output 'dvect' has the length of N-nsep, where N is length of the input
% 'vect', but note some
% elements in 'sumvect' does not have the same number of 'nsum' contributing
% to the sum, depending on the 'direction'.
%
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/06/09
% Last modified date:   2022/06/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('direction','forward');

% vect = reshape(vect,[],1);
N = size(vect,1);
dvect = zeros(N-nsep, size(vect,2));

if isequal(direction,'forward')
  for i = 1+nsep: N 
    dvect(i-nsep, :) = vect(i, :)-vect(i-nsep, :);
  end
elseif isequal(direction,'backward')
  for i = 1: N-nsep 
    dvect(i, :) = vect(i+nsep, :)-vect(i, :);
  end
end