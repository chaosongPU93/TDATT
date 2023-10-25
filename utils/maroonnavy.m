function c = maroonnavy(m)
%MAROONNAVY    Shades of maroon red and navy blue color map
%   MAROONNAVY(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(marronnavy)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

% Chao Song, chaosong@princeton.edu
% Inspired by REDBLUE.m 
% First created date:   2021/10/21
% Last modified date:   2021/10/21

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 0.5] to [1 1 1], then [1 1 1] to [0.5 0 0];
    m1 = m*0.5;
    r1 = (0:m1-1)'/max(m1-1,1);
    g1 = r1;
%     if (mod(m1,2) == 0)
%       m2 = m1*0.5;
      b1 = (m1-1:m-2)'/max(m-2,1);
%       b1 = (m2-1:m1-2)'/max(m1-2,1);
      r2 = flipud(b1);
      g2 = flipud(r1);
      b2 = flipud(g1);
      r = [r1; r2];
      g = [g1; g2];
      b = [b1; b2];
%     else
%       m2 = ceil(m1*0.5);
%       b1 = (m2:m1)'/max(m1+1,1);
%       r2 = flipud(b1);
%       g2 = flipud(r1);
%       b2 = flipud(g1);
%       r = [r1; 1; r2];
%       g = [g1; g2];
%       b = [b1; 1; b2];
% 
%     end
else
    % From [0 0 0.5] to [1 1 1] to [0.5 0 0];
    m1 = floor(m*0.5);
    r1 = (0:m1-1)'/max(m1,1);
    g1 = r1;
%     if (mod(m1,2) == 0)
%       m2 = m1*0.5;
      b1 = (m1:2*m1-1)'/max(2*m1,1);
%       b1 = (m2-1:m1)'/max(m1,1);
      r2 = flipud(b1);
      g2 = flipud(r1);
      b2 = flipud(g1);
      r = [r1; 1; r2];
      g = [g1; 1; g2];
      b = [b1; 1; b2];
%     else
%       m2 = ceil(m1*0.5);
%       b1 = (m2:m1)'/max(m1+1,1);
%       r2 = flipud(b1);
%       g2 = flipud(r1);
%       b2 = flipud(g1);
%       r = [r1; 1; 1; r2];
%       g = [g1; 1; 1; g2];
%       b = [b1; 1; 1; b2];
%     end
end

c = [r g b]; 

