% function compare_bursts_at_diffreg(FLAG,TYPE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function compare_bursts_at_diffreg(FLAG,TYPE)
% -- This script is to plot the comparison of tremor bursts active in each 
% of the ETS episode at different regions to their relationship, e.g.,
% overlapping in time; length in time; dominancy.
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/10/26
% Last modified date:   2021/10/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


format short e   % Set the format to 5-digit floating point
clc
clear
close all

%% default value for easy debugging
defval('FLAG','PGC');
defval('TYPE','HF');

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

% get the scrsz in pixels and number of pixels per inch of monitor 1
[scrsz, res] = pixelperinch(1);

workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
hypopath = strcat(getenv('MHOME'), '/Seisbasics/hypoinverse/forproj21');

% FLAG = 'PGC';
% FLAG = 'TWKB';

disp(FLAG);

if strcmp(FLAG, 'TWKB')     % use TWKB catalog
  temppath = strcat(datapath, '/templates/LZBtrio/');
  rstpath = strcat(datapath, '/LZBtrio');
  figpath = strcat(workpath, '/project2021/TWKBtrio');
  
  nfampool = [
              '002';
              '043';
              '141';
              '047';
              '010';
              '144';
              '099';
              '068';
              '125';
              '147';
              '017';
              '006';
              '001';
              ];
  nfam = size(nfampool,1);
  
  stas=['TWKB '
    'LZB  '
    'MGCB '];     % determine the trio and order
  nsta=size(stas,1);         %  number of stations
  
  
elseif strcmp(FLAG, 'PGC')     % use PGC catalog
  temppath = strcat(datapath, '/templates/PGCtrio/');
  rstpath = strcat(datapath, '/PGCtrio');
  figpath = strcat(workpath, '/project2021/PGCtrio');
  
  nfampool = [
                '002';
                '243';
                '240';  % the most recent catalog of PGC is from fam 002, 243 and 240
%                 '253';
%                 '036';
%                 '251';
              ];
  nfam = size(nfampool,1);
  
  stas=['PGC  '
    'SSIB '
    'SILB '
    'KLNB '
    ];
  nsta=size(stas,1);         %  number of stations
  
end


%%% load the time range of the tremor bursts
%%% option 1, windows from automatic clustering algorithm
bndflag = ['bnd   '
           'bndse '
           'bndcnt'
           'bndnw '
           'bndn  '
           'bndne '
           ];
nreg = size(bndflag,1);
for i = 1: nreg
  fname = strcat(hypopath, '/tremor_burst_ranges_',FLAG,TYPE,strtrim(bndflag(i,:)));
  temp = load(fname);
  dates = unique(temp(:,1));
  years = unique(floor(dates/1000));
  nets = length(years);
  for j = 1: nets
    year = years(j);
    temp2 = temp(floor(temp(:,1)/1000)==year, :);
    burstcell{i,j} = temp2;
  end
end

% %%% option 2, windows from Allan's eyes using ginput
% jdays= [2003 062; %jdays are for data.
%         2003 063;
%         2004 196;
%         2004 197;
%         2005 254;
%         ];
% ntrange = [];      
% for id = 1: size(jdays,1)        
%   wins=getspecwins(jdays(id, :));
%   ran = zeros(size(wins,1), 3);
%   ran(:, 1) = (jdays(id, 1)*1000 + jdays(id, 2))*ones(size(wins,1),1);
%   ran(:, 2:3) = wins;
%   ntrange = [ntrange; ran];
% end

% keyboard

%% plot the comparison of tremor bursts active at different regions, for different ETSs
for j = 1: nets
  for i = 1: nreg
%     %%% use all bursts    
%     burstets{i,1} = burstcell{i,j};

    %%% use only bursts that are longer than the length needed by spectrum computation
    temp = burstcell{i,j};
    temp2 = [];
    count = 0;
    for k = 1: size(temp,1)
      tst = temp(k,2);
      ted = temp(k,3);
      wlen = 2048;
      sps = 100;
      wlensec = wlen/sps;
      if ted-tst>=wlensec
        count = count+1;
        temp2(count,:) = temp(k,:);
      end
    end
    burstets{i,1} = temp2;

  end
  
  % plot bursts at different regions
  [f] = plt_burst_diffreg(50,burstets);
  for i = 1: nreg 
    text(f.ax(i),0.04,0.9,strtrim(bndflag(i,:)),'FontSize',11,'unit','normalized',...
      'horizontalalignment','left','EdgeColor','k','Margin',2);
  end
  print(f.fig,'-dpdf',strcat(figpath,'/burstdiffreg_',FLAG,TYPE,num2str(years(j)),'.pdf'));

end


