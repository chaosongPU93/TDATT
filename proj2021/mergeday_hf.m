function mergeday_hf(fam,fcflag,convflag,iup)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the funcion to combine the detections of all days in any fam
% using TWKB trio, similar to 'mergeday_hf_PGC'
%
%
% NOTE:
%   filtering effect correction should be
%   done in this step, i.e., set the fcflag to be 1
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/07/24
% Last modified date:   2021/07/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('fam', '043');
defval('fcflag', 1);
defval('convflag', 0);
defval('iup', 4);

%% Initialization
%%% SAME if focusing on the same region (i.e. same PERMROTS and POLROTS)
%%% AND if using the same family, same station trio

format short e   % Set the format to 5-digit floating point
% clear
% close all

% set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

% WHEN CHANGING FAMILIES CHANGE: (1)dates (2)Bostnames (3)hilo,frequency band 
% (4)mshift (5)bostsec (6)stas (7)PERMROTS and POLROTS (8)tempoffs

workpath = getenv('ALLAN');
%%% choose the data path here 
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
rstpath = strcat(datapath, '/LZBtrio');

% fam='243';     % family number
freqflag='hf';  % flag to indicate whether to do hf or lf;

FLAG = 'TWKB';

timoffrot= [
            2003 060;
            2003 061;
            2003 062; % dates of the year, two columns: year Julian_day
            2003 063;
            2003 064;   % LZB only has data from 2003 060-064
            2004 194;
            2004 195;
            2004 196;
            2004 197;
            2004 198;
            2004 199;
            2004 200;
            2004 201;
            2004 202;
            2004 203;
            2005 254;
            2005 255;
            2005 256;
            2005 257;
            2005 258;
            2005 259;
            2005 260;
            2005 261];

nday = size(timoffrot, 1);

%%% corresponds to PERMROTS
PERMSTA=['PGC'        % permanent station names
         'LZB'];
POLSTA=['SSIB '           % polaris station names
        'SILB '
        'KLNB '
        'MGCB '
        'TWKB '];
    
stas=['TWKB '
      'LZB  '
      'MGCB '];     % determine the trio and order
  
%% Important parameters same as that during detections    
nsta=size(stas,1);         %  number of stations
sps=40;     % samples per second

%%% IMPORTANT, NEED to change sometimes
%Basics of the cross-correlation:  Window length, number of windows, filter parameters, etc.
winlensec=4;
winoffsec=1;        % window offset in sec, which is the step of a moving window
winlen=winlensec*sps;      % length in smaples
hi=6.5;    % frequency band
lo=1.25;
npo=2;     % poles, passes of filters
npa=2;
cyclskip = 0;
mshift=14+cyclskip; %19; %maximum shift for the x-correlations. 19 for 002 Stanford,    % in sps, 0.5s*40sps=20
loopoffmax=1.5; %1.5 for standard 1.5-6Hz; 4 for 0.5-1.5Hz.  2 for non-interpolated.   % what is loopoffmax, the circuit of time offsets
xcmaxAVEnmin=0.4; %0.44; %0.44 for 002 Stanford %0.45; %0.36 for 4s 1-12 Hz; %0.4 for 4s 1.5-6 Hz and 6s 0.5-1.5Hz; 0.36 for 4s2-8 Hz ; 0.38 for 4s0.75-6 Hz; 0.094 for 128s 2-8Hz;  0.1 for 128s 1.5-6Hz; 0.44 for 3-s window?
    
if fcflag
    %%% load the template filtering effect result
    lolf = 0.5;
    hilf = 1.25;
    lohf = 1.25;
    hihf = 6.5;
    winsechf = winlensec;
    winseclf = 16;
    lofflf = 4;
    ccminlf = 0.35;
    PREFIX = strcat(fam,'.lf',num2str(lolf),'-',num2str(hilf),'.hf',num2str(lohf),'-',num2str(hihf),...
        '.hw',num2str(winsechf),'.lw',num2str(winseclf),'.lo',num2str(lofflf),'.ccm', ...
        num2str(ccminlf),'.','80sps');
    fname = strcat(rstpath, '/MAPS/tempfeff_LZB_enc','_',PREFIX);
    filtcor = load(fname);
    spsratio = sps/80;
    filthf = filtcor(:,1)*spsratio;   % lf template shift due to filterig effect, sign is the same
else
    filthf = zeros(2,1);
end
  
  
%% START TO LOOP FOR every day
rstall = [];
%cycle over each day:
for nd=1:length(timoffrot(:,1))      % num of rows, also num of days
%     nd=1;

    year=timoffrot(nd,1);
    YEAR=int2str(year);
    jday=timoffrot(nd,2);
    if jday <= 9
        JDAY=['00',int2str(jday)];
    elseif jday<= 99
        JDAY=['0',int2str(jday)];
    else
        JDAY=int2str(jday);
    end
    
    IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),...
            '.',int2str(npo),int2str(npa),'.ms', int2str(mshift)]
        
    if iup == 1     
        fname1 = strcat(rstpath,'/MAPS/pj21mapall',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '4add');   
        fname2 = strcat(rstpath,'/MAPS/pj21mapall',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '3add');
    else
        fname1 = strcat(rstpath,'/MAPS/pj21mapallup',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '4add');   
        fname2 = strcat(rstpath,'/MAPS/pj21mapallup',IDENTIF,'_',num2str(lo),'-',num2str(hi),'_',...
                        num2str(winlen/sps),'s',num2str(sps),'sps', '3add');
    end
    
    if isfile(fname1)
        rstday = load(fname1);
        %%% structure of allrst: 30+4*nstanew cols, if 4 new stas, then will be 46 cols
        %%% UPDATED at 2021/03/20
        %%%   1:timswin(n) 2:xmaxSTA12ntmp(n) 3:xmaxSTA13ntmp(n) 4:xcmaxAVEnbang(nin) 5:loopoff(n)
        %%%   6:cumsumtrdiff 7:timswin(n)-winlensec/2+idiff/sps 8:cumsumtrdiff/cumsumtr(winlen)
        %%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
        %%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
        %%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
        %%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
        %%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
        %%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)

        %%% additional stations for checking detections
        stasnew=['PGC  '
                 'SSIB '
                 'SILB '
                 'KLNB '];  % twkb lzb mgcb   
        if year==2003   % in 2003, station KLNB is named KELB, so use KELB to replace KLNB in 2003
            POLSTA(3,:)='KELB ';
            stasnew(4,:)='KELB ';
        else
            POLSTA(3,:)='KLNB ';  % remember to change it back
            stasnew(4,:)='KLNB ';
        end
        nstanew=size(stasnew,1);
    elseif isfile(fname2)
        rstday = load(fname2);
        stasnew=['PGC  '
                 'SSIB '
                 'SILB '];
        nstanew=size(stasnew,1);
    else
        fprintf('Day %s %s of fam %s is skipped because of no detections.\n',YEAR, JDAY, fam);
    end   
  
    %%% filtering correction
    if size(filthf,1) > 2   % which means containing the effects of additional stations
        rstday(:, 2) = rstday(:, 2) - filthf(2);  % the original correction contains all seven with diff algorithm
        rstday(:, 3) = rstday(:, 3) - filthf(3);  % ranges with order
    else
        rstday(:, 2) = rstday(:, 2) - filthf(1);  % since the total offset 12 contains filtering effect, so substract to correct
        rstday(:, 3) = rstday(:, 3) - filthf(2);
    end
    
    %%% whether to convert offsets from different fams to a uniform frame
    %%% ALGORITHM is confirmed to be true!
    CATA = 'new';
    [PERMROTS, POLROTS] = GetRotsCommon(FLAG,fam,CATA,datapath,0,0);
    reftfam = POLROTS(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
    PERMROTS(:,4) = PERMROTS(:,4)-reftfam;    % to make sure that 1st station is 0
    POLROTS(:,4) = POLROTS(:,4)-reftfam;

    reffam = '043';     %043
    [refperm, refpol] = GetRotsCommon(FLAG,reffam,CATA,datapath,0,0);
    reftref = refpol(5,4);     % reftime is the reference time at the 1st station,depends on the choice of 1st station
    refperm(:,4) = refperm(:,4)-reftref;    % to make sure that 1st station is 0
    refpol(:,4) = refpol(:,4)-reftref;
  
    % conversion time constants, NOTE that the slow/fast time offset and offset between stas
    % is all calculated based on 40 sps
    conv(1) = PERMROTS(2,4) - refperm(2,4);   % station 2, LZB
    conv(2) = POLROTS(4,4) - refpol(4,4);    % station 3, MGCB
    if convflag
        disp('Offsets from different fams are converted to a uniform frame ...');
        spsrat = sps/40;
        rstday(:, 2) = rstday(:, 2)-conv(1)*spsrat;   % get the new off based on the ref fam frame
        rstday(:, 3) = rstday(:, 3)-conv(2)*spsrat;
    end
    
    STA12offsec = rstday(:, 2)/sps;
    STA13offsec = rstday(:, 3)/sps;

    famarr = ceil(str2double(fam))*ones(size(rstday,1), 1); % label of fam number array
    datearr = (year*1000+jday)*ones(size(rstday,1), 1);

    %reorder the columns and so on to make the structure be as follows:
    %%% 34+4*nstanew cols, if 4 new stas, then will be 50 cols
    %%% UPDATED at 2021/06/23
    %%%   1:off12(n) 2:off13(n) 3:off12sec 4:off13sec 5:fam 6:date 7:timswin(n)
    %%%   8:timswin(n)-winlensec/2+idiff/sps 9:xcmaxAVEnbang(nin) 10:loopoff(n)
    %%%   11:cumsumtrdiff  12:cumsumtrdiff/cumsumtr(winlen)
    %%%   for the rest, n=n+4;
    %%%   9:Ampsq(nin) 10:Prior(nin) 11:Post(nin) 12:Prior4(nin) 13:Post4(nin) 14:Prior8(nin)
    %%%   15:Post8(nin) 16:Prior16(nin) 17:Post16(nin) 18:cc(nin) 19:ccprior(nin)
    %%%   20:ccpost(nin) 21:ccprior4(nin) 22:ccpost4(nin) 23:STAamp(nin,2) 24:STAamp(nin,3)
    %%%   25:sumsSTA12n(n,iSTA12bang) 26:sumsSTA13n(n,iSTA13bang)
    %%%   27:sumsSTA32n(n,iSTA32bang) 28:sigma(nin) 29:sigma12(nin) 30:sigma13(nin)
    %%%   in(1:nstanew) loff(1:nstanew) ioff(1:nstanew) ccmaxave(1:nstanew)

    rstsave = [rstday(:, 2:3) STA12offsec STA13offsec famarr datearr rstday(:, 1) rstday(:, 7) ...
               rstday(:, 4:6) rstday(:, 8:end)];

    % concatenate them all
    rstall = [rstall; rstsave];
  
end  


%%% save to file
if convflag == 0 
    if iup == 1
        PREFIX = strcat(fam,'.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',int2str(npo),...
                        int2str(npa),'.ms', int2str(mshift));
    else
        PREFIX = strcat(fam,'.up.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',int2str(npo),...
                        int2str(npa),'.ms', int2str(mshift));
    end                
else
    if iup == 1
        PREFIX = strcat(fam,'.conv.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));
    else
        PREFIX = strcat(fam,'.upconv.lo',num2str(loopoffmax),'.cc',num2str(xcmaxAVEnmin),'.',...
                        int2str(npo),int2str(npa),'.ms', int2str(mshift));
    end
end

fid = fopen(strcat(rstpath, '/MAPS/pj21timeori_',PREFIX,'_',num2str(lo),'-',num2str(hi),'_',...
            num2str(winlen/sps),'s',num2str(sps),'sps', num2str(nstanew), 'add'),'w+');
if nstanew == 3
    fprintf(fid,'%6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %.4f %.4f %.4f %d %d %d %.3f %.3f %.3f \n',...
        rstall');
elseif nstanew == 4
    fprintf(fid,'%6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        rstall');
end
fclose(fid);

% keyboard       

       
       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
  
  