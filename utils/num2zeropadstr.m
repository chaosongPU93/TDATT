function str = num2zeropadstr(num,len)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% str = num2zeropadstr(num,len)
%
% Sometimes we want transform a bunch of number to strings which have the same
% length to easier processing. This function just uses 'num2str' to turn a
% number to string first, then pad zeros to the start, to make it has a fixed
% length.
%
% 
% By Chao Song, chaosong@princeton.edu
% First created date:   2022/04/04
% Last modified date:   2022/04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = reshape(num,[],1);
nnum = size(num,1);
% str = string(nnum, len);

for i = 1: nnum
  if len == 2
    if num(i) < 10
      str(i,1:len) = strcat(repmat('0',1,len-1),num2str(num(i)));
    else
      str(i,1:len) = num2str(num(i));
    end
    
  elseif len == 3
    if num(i) < 10
      str(i,1:len) = strcat(repmat('0',1,len-1),num2str(num(i)));
    elseif num(i) < 1e2
      str(i,1:len) = strcat(repmat('0',1,len-2),num2str(num(i)));
    else
      str(i,1:len) = num2str(num(i));
    end
    
  elseif len == 4
    if num(i) < 10
      str(i,1:len) = strcat(repmat('0',1,len-1),num2str(num(i)));
    elseif num(i) < 1e2
      str(i,1:len) = strcat(repmat('0',1,len-2),num2str(num(i)));
    elseif num(i) < 1e3
      str(i,1:len) = strcat(repmat('0',1,len-3),num2str(num(i)));
    else
      str(i,1:len) = num2str(num(i));
    end
    
  elseif len == 5
    if num(i) < 10
      str(i,1:len) = strcat(repmat('0',1,len-1),num2str(num(i)));
    elseif num(i) < 1e2
      str(i,1:len) = strcat(repmat('0',1,len-2),num2str(num(i)));
    elseif num(i) < 1e3
      str(i,1:len) = strcat(repmat('0',1,len-3),num2str(num(i)));
    elseif num(i) < 1e4
      str(i,1:len) = strcat(repmat('0',1,len-4),num2str(num(i)));
    else
      str(i,1:len) = num2str(num(i));
    end
    
  elseif len == 6
    if num(i) < 10
      str(i,1:len) = strcat(repmat('0',1,len-1),num2str(num(i)));
    elseif num(i) < 1e2
      str(i,1:len) = strcat(repmat('0',1,len-2),num2str(num(i)));
    elseif num(i) < 1e3
      str(i,1:len) = strcat(repmat('0',1,len-3),num2str(num(i)));
    elseif num(i) < 1e4
      str(i,1:len) = strcat(repmat('0',1,len-4),num2str(num(i)));
    elseif num(i) < 1e5
      str(i,1:len) = strcat(repmat('0',1,len-5),num2str(num(i)));
    else
      str(i,1:len) = num2str(num(i));
    end
    
  else
    error('Length of zero-padded string is longer than the max allowed by the function \n');
  end
  
end  
  
% keyboard  
  
  
  
  