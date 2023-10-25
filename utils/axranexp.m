function varargout=axranexp(ax,xory,perc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax = axranexp(ax,xory,perc)
%
% This function is to expand or shrink the range of the axis by some 
% percentage. 'xory' determine which limit is going to be changed. 
% Similar in usage to FJ Simons' function 'openup.m'
%
% INPUT:
% 
%  ax       Axis handles (default: gca)
%  xory     1 right
%           2 top [default]
%           3 left
%           4 bottom
%           5 left and right
%           6 top and bottom
%  perc     Percentage of the data range [default: 10], positive means
%           expanding, negative means shrinking
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/03/12
% Last modified date:   2022/03/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('ah',gca);
defval('xory',2);
defval('perc',10);

switch xory
 case {1,3,5}
  wat='xlim';
 case {2,4,6}
  wat='ylim';
end

for i = 1: length(ax) 
  
  switch xory
    case {1,3,5}
      xran = ax(i).XLim;
      xlen = xran(2)-xran(1);
      if xory == 1
        xmax = xlen*(1+perc/100)+xran(1);
        nxran = [xran(1) xmax];
      elseif xory == 3
        xmin = xran(2)-xlen*(1+perc/100);
        nxran = [xmin xran(2)];
      else
        xmin = xran(2)-xlen*(1+perc/100);
        xmax = xlen*(1+perc/100)+xran(1);
        nxran = [xmin xmax];
      end
      set(ax(i),wat,nxran);
    case {2,4,6}
      yran = ax(i).YLim;
      ylen = yran(2)-yran(1);
      if xory == 2
        ymax = ylen*(1+perc/100)+yran(1);
        nyran = [yran(1) ymax];
      elseif xory == 4
        ymin = yran(2)-ylen*(1+perc/100);
        nyran = [ymin yran(2)];
      else
        ymin = yran(2)-ylen*(1+perc/100);
        ymax = ylen*(1+perc/100)+yran(1);
        nyran = [ymin ymax];
      end
      set(ax(i),wat,nyran);     
  
  end
  
end

% Optional output
varns={ax};
varargout=varns(1:nargout);






