function newd = kdensify1d(oldd, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the script to linearly densify the 1d vector k times denser, say
% the original vector size is (n,1), after this densification, the size would
% be (k(n-1)+1, 1)
% 
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/04/08
% Last modified date:   2020/04/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k=4;
oldd = reshape(oldd,[],1);
odelta = oldd(2:end)-oldd(1:end-1);
ndelta = odelta/4;
newd = [];
for i = 1: size(oldd,1)-1
    tmp = linspace(oldd(i),oldd(i+1)-ndelta(i), k);
    newd = [newd; tmp'];
end
newd = [newd; oldd(end)];