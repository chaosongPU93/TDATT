% function [respfiles] = fetchRESP(filelist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is a function to download the RESP files from IRIS DMC by creating a 
% query url and done by 'wget' command. 
%
%
%   INPUT:
%       f:      the construct containing frame and for figure
%       d1:     the start day of the time axis, Must be datestring format
%       d2:     the end day of the time axis, Must be datestring format
%       ievtday:    the regular event day on the time axis, counting from the start 
%       nevtobj:    the event counts on that day
%
%   OUTPUT:
%       f:    the construct for figure
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/02/20
% Last modified date:   2020/02/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% system('')

fid = fopen('/home/data2/chaosong/shikoku_kii/20141004/filelist');
% Read all lines & collect in cell array
txt = textscan(fid,'%s','delimiter','\n');

ntrace = size(txt{1});
for i = 1: ntrace
    trnm = char(txt{1}(i));
    trnmsp = strsplit(trace,'.');
    
    ntwk = char(trnmsp{1});
    stnm = char(trnmsp{2});
    comp = char(trnmsp{3});
end