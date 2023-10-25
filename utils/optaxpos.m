function varargout = optaxpos(f,nrow,ncol,xran,yran,xsep,ysep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axpos = optaxpos(f,nrow,ncol,xran,yran,xsep,ysep)
%
% This function is to return an optimal axis positions of subplots with an 
% equal length and width, based on the input of number of rows and cols, 
% range for plotting in x and y directions, separation of subplots in 
% x and y directions, in the normlized length units (a total of 1)
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/29
% Last modified date:   2021/11/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('xran',[0.1 0.9]);
defval('yran',[0.1 0.9]);
defval('xsep',0.05);
defval('ysep',0.05);

axpos = zeros(nrow*ncol,4);
for irow = 1: nrow
  for icol = 1: ncol
    isub = (irow-1)*ncol+icol;
    xwid = (xran(2)-(ncol-1)*xsep-xran(1))/ncol;
    ywid = (yran(2)-(nrow-1)*ysep-yran(1))/nrow;
    xloc = (icol-1)*(xsep+xwid)+xran(1);
    yloc = (nrow-irow)*(ysep+ywid)+yran(1);
    axpos(isub,:) = [xloc yloc xwid ywid];
  end
end

for isub = 1:nrow*ncol
  set(f.ax(isub), 'position', axpos(isub,:));
end

% Optional output
varns={axpos};
varargout=varns(1:nargout);