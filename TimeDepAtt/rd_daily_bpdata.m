function [STAopt,STAort,STAnzeros,fileflag] = rd_daily_bpdata(year,jday,prename,stas,PERMSTA,POLSTA,...
          PERMROTS,POLROTS,sps,lo,hi,npo,npa,nwin,winlen,winoff,igstart,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to read the daily horizontal seismograms with shear-wave splitting
% corrected, particle motion oriented, and bandpassed with a specified
% frequency range. Allows the separation of a full day into sub-windows
% with overlapping etc.
% Although it takes in the absolute path of the data file, it has the 'factor'
% hard-coded which is used to manually account for the gain difference between
% the PGC and LZB, and between Polaris stations before and after 2003213,
% to make sure they would have a similar amplitude level.
%
% --2022/11/30, I changed the 'fact' at KLNB to be 1.5e-3 for ALL time. Note
%   that the seismogram is only scaled by 'fact', no 'scaleseisms' used
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/11/08
% Last modified date:   2021/11/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

YEAR = int2str(year);
if jday <= 9
  JDAY=['00',int2str(jday)];
elseif jday<= 99
  JDAY=['0',int2str(jday)];
else
  JDAY=int2str(jday);
end
MO=day2month(jday,year);     % EXTERNAL function, day2month, get the month of one particular date
% %   IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
% %     int2str(npo),int2str(npa),'.ms', int2str(mshift)]
% direc=[path, '/arch', YEAR,'/',MO,'/'];     % directory name
% prename=[direc,YEAR,'.',JDAY,'.00.00.00.0000.CN'];    %  path plus prefix of data file,

if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
  POLSTA(3,:)='KELB ';
else
  POLSTA(3,:)='KLNB ';  % remember to change it back
end

%in case you mix-use the 'KLNB' or 'KELB' in the station name, this would correct it for you
if ismember('KLNB ',stas,'rows') || ismember('KELB ',stas,'rows')
  [~,idx1]=ismember('KELB ',stas,'rows');
  [~,idx2]=ismember('KLNB ',stas,'rows');
  if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
    stas(idx1+idx2,:)='KELB ';
  else
    stas(idx1+idx2,:)='KLNB ';  % remember to change it back
  end
end


%Read the data; find glitches (many consecutive zeros) and flag those windows in STAnzeros.
%Get timsSTA from the permanent stations (last one over-writes):
tracelen=86400*sps; %one day of data at 40 sps, overall trace length, 24*3600
nsta = size(stas,1);
STAopt=zeros(tracelen,nsta+1);    % one day of samples
STAort=zeros(tracelen,nsta+1);
STAnzeros=zeros(nwin,nsta);    % used to flag windows with glitches

fileflag = 1;   % 1 means all files exist, the file status is normal
for ista=1:nsta
  found=0;
  fname = [];
  
  [LIA,idx]=ismember(stas(ista,:),PERMSTA,'rows');     % ismember, to determine whether each row of stas is contained in PERMSTA, return logical value 1/0 and index
  if LIA      % if station is a permanent one
    found=found+LIA;
    if strcmp(PERMSTA(idx,1:3),'PGC')     % string compare
      %%% WAHT is the meanning of 'fact', similar to instrument
      %%% response
      fact=1.0e-3;
    elseif strcmp(PERMSTA(idx,1:3),'LZB')
      fact=1.8e-3;
    end
    %%% readperms is an EXTERNAL FUNCTION
    % opt: optimal seismogram after rotations
    % ort: orthogonal seismogram after rotations
    % nzeros: number of zeros in the trace
    % timsSTA: time sequence
    fname = strcat(prename,'.',PERMSTA(idx,1:3),'..BHE.D.SAC');
    if isfile(fname)    % if have the data file
      %                 % 1. this is for data without removing station response
      %                 [opt,ort,nzeros,timsSTA]=readperms(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
      % 2. this is for data with no response
      [opt,ort,nzeros,time]=readpermsv2(prename,PERMSTA,PERMROTS,idx,sps,lo,hi,npo,...
        npa,fact,nwin,winlen,winoff,igstart,varargin);
      STAopt(:, 1)=time;
      STAort(:, 1)=time;  
    else
      fileflag = 0;   % change the file flag to 0, meaning abnormal
      fprintf('No data for station %s in day %s %s, this day will be omitted. \n',...
        PERMSTA(idx,1:3), YEAR, JDAY);
      break   % break the entire station loop
    end
  end
  
  [LIA,idx]=ismember(stas(ista,:),POLSTA,'rows');
  if LIA        % if are in POLSTA
    found=found+LIA; %better be 1
    if year==2003 && jday<213 && ~strcmp(POLSTA(idx,1:4),'KELB')      % should result from some criteria
%     if year==2003 && jday<213     % should result from some criteria
      fact=7.5e-3;
    elseif year==2003 && jday<213 && strcmp(POLSTA(idx,1:4),'KELB')      % should result from some criteria
      fact=5.0e-3;
    else
      fact=1.5e-3;
    end
    
    fname = strcat(prename,'.',POLSTA(idx,1:4),'..HHE.D.SAC');
    if isfile(fname)    % if have the data file
      %                 % 1. this is for data without removing station response
      %                 [opt,ort,nzeros]=readpols(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,nwin,winlen,winoff,igstart);
      % 2. this is for data with no response
      [opt,ort,nzeros]=readpolsv2(prename,POLSTA,POLROTS,idx,sps,lo,hi,npo,npa,fact,...
        nwin,winlen,winoff,igstart,varargin);
    else
      fileflag = 0;   % change the file flag to 0, meaning abnormal
      fprintf('No data for station %s in day %s / %s. \n',POLSTA(idx,1:4), YEAR, JDAY);
      break   % break the entire station loop
    end
  end

  if fileflag
    STAopt(:, ista+1)=opt;
    STAort(:, ista+1)=ort;
    STAnzeros(:, ista) = nzeros;
  end
end


% keyboard


