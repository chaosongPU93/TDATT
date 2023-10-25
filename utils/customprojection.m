function [projx,projy,projxy] = customprojection(locxy,projang)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [propx,orty,nlocxy] = customprojection(olocxy,propang)
%
% Sometimes you might have the need to project the location coordinate
% of sources to along a custom direction. Basically it is a coordinate 
% frame transformation, ie., get the location in a custom frame. Call
% the 'coordinate_rot' for actual transform.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2023/01/26
% Last modified date:   2023/01/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projxy = locxy;
for i = 1: size(locxy,1)
  [projxy(i,1), projxy(i,2)] = coordinate_rot(locxy(i,1),locxy(i,2),-(projang-90),[0 0]);
end
projx = projxy(:,1);
projy = projxy(:,2);