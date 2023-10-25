%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is the script to analyse the inversion results from hypoinverse, set
% a cutoff distance, e.g. 8 km. from the center of lfe fam to their own 
% detections. Perserve the ones within this distance only.
%
%   version for fams using TWKB trio, currently for our purpose, it is not
%   necessary that the distance cutoff has to be done. So at present, this
%   is served as analyzing the inversion results mainly
%
%
% Modified by Chao Song, chaosong@princeton.edu
% First created date:   2021/07/24
% Last modified date:   2021/07/24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
format short e   % Set the format to 5-digit floating point
clear
close all

set(0,'DefaultFigureVisible','on');
% set(0,'DefaultFigureVisible','off');   % switch to show the plots or not

scrsz=get(0,'ScreenSize');

workpath = getenv('MHOME');
rstpath = strcat(workpath, '/Seisbasics/hypoinverse/forproj21');

winlenhf = 4;

winlenlf = 16;
lofflf = 4;
ccminlf = 0.35;
        
        
%% for detections in time
SUFFIXhf = strcat('up.hf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.lzbfam.pj21.',SUFFIXhf);
hfmaptime = load(fname);
% 53 cols, format is:
%   lon lat dep + previous 50 cols
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



SUFFIXlf = strcat('lf.time.',num2str(winlenhf),'_',num2str(winlenlf),'.',num2str(lofflf),...
                  '.',num2str(ccminlf));
fname = strcat(rstpath, '/evtloc.lzbfam.pj21.',SUFFIXlf);
lfmaptime = load(fname);

%%% for the locations of lfe families, by inversion from hypoinverse
%%% this is inverted from (0,0) of all fams, same order, location of control points
%%% inverted from /home/data2/chaosong/Seisbasics/hypoinverse/lzbrst/chat00p.reffam
loccont = [-123.492667 48.451500 38.1400; 
           -123.772167 48.493000 35.5900; 
           -123.863167 48.528167 35.2100;
           -123.603333 48.440167 36.7100;
           -123.800167 48.408833 34.5200;
           -123.893333 48.536500 35.0700;
           -123.864500 48.498667 34.8800;
           -123.753333 48.525667 36.2000;
           -123.703667 48.502667 36.4100;
           -123.814333 48.538667 35.7900;
           -123.838500 48.544833 35.6600;
           -123.908000 48.494167 34.5100;       % 006
           -123.879667 48.446167 34.2600;       % 001
%            -123.967000 48.458667 33.7800;       % 158, 20200916,testing purpose
           ];

relacont = loccont;
% [dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),loccont(2,1),loccont(2,2));
[dx, dy] = absloc2relaloc(loccont(:,1),loccont(:,2),-123.772167,48.493000);
relacont(:,1) = dx;
relacont(:,2) = dy;       

% convert absolute loc to relative loc to its own lfe fam
nfampool = ['002';
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
disp(nfam);


% %% read the alligned seismograms of all detections and get the amplitude info
% %%% NOTE: it is impossible to track the order of detections since they have been sorted somewhere
% %%% thus the only way to link them is the fam number, day and time
% timoffrot= [
%             2003 061
%             2003 062;
%             2003 063;
%             2004 196;
%             2004 197;
%             2004 198;
%             2004 199;
%             2005 254;
%             2005 255;
%             2005 256];
% nday = size(timoffrot, 1);
% 
% tracepath = '/home/data2/chaosong/matlab/allan/data-no-resp/PGCtrio';
% amphf = zeros(size(hfmaptime,1),2);   % hf, amplitude of detections for all days of all fams
% for i=1: size(nfampool,1)
%     fam = nfampool(i,:);
%     for nd=1:length(timoffrot(:,1))      % num of rows, also num of days
%         year=timoffrot(nd,1);
%         YEAR=int2str(year);
%         jday=timoffrot(nd,2);
%         if jday <= 9
%             JDAY=['00',int2str(jday)];
%         elseif jday<= 99
%             JDAY=['0',int2str(jday)];
%         else
%             JDAY=int2str(jday);
%         end
%         IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo2.1.cc0.44.22.ms29_1.25-6.5_4s40sps4add'];
%         fname = strcat(tracepath, '/MAPS/pj21traceallup_',IDENTIF);
%         if isfile(fname)
%             dettrace = load(fname);
%             ndet = round(size(dettrace,1)/(4*40));
%             for j = 1: ndet
%                 ampmax = mean(max(dettrace((j-1)*4*40+1:j*4*40, 2:4)));
%                 ampmin = mean(min(dettrace((j-1)*4*40+1:j*4*40, 2:4)));
%                 midtime = dettrace((j-1)*4*40 +1 +4*40/2, 1);
%                 midtime = sprintf('%.1f',midtime);
%                 midtime = str2double(midtime);
%                 objind = find(hfmaptime(:,8)==str2double(fam) & ...
%                               hfmaptime(:,9)==str2double(strcat(YEAR,JDAY)) & ...
%                               abs(hfmaptime(:,10)-midtime)< 0.01);
%                 if isempty(objind) || size(objind,1) > 1
%                     disp(midtime)
%                     disp(strcat(YEAR,JDAY))
%                     disp(objind)
%                     break
%                 end
%                 amphf(objind,1) = ampmax; 
%                 amphf(objind,2) = ampmin;
%             end
%         else
%             continue
%         end
%     end
% end
% disp('Amplitude linked to hf detections done');
% 
% 
% amplf = zeros(size(lfmaptime,1),2);   % lf, amplitude of detections for all days of all fams
% for i=1: size(nfampool,1)
%     fam = nfampool(i,:);
%     for nd=1:length(timoffrot(:,1))      % num of rows, also num of days
% % nd=2;
%         year=timoffrot(nd,1);
%         YEAR=int2str(year);
%         jday=timoffrot(nd,2);
%         if jday <= 9
%             JDAY=['00',int2str(jday)];
%         elseif jday<= 99
%             JDAY=['0',int2str(jday)];
%         else
%             JDAY=int2str(jday);
%         end
%         if isequal(fam,'002')
%             IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo4.cc0.5.22.ms39_0.5-1.25_16s20sps4add'];
%         else
%             IDENTIF=[YEAR,'.',JDAY,'.',fam,'.lo4.cc0.45.22.ms39_0.5-1.25_16s20sps4add'];
%         end
%         fname = strcat(tracepath, '/MAPS/pj21traceall_',IDENTIF);
%         if isfile(fname)
%             dettrace = load(fname);
%             ndet = round(size(dettrace,1)/(16*20));
%             for j = 1: ndet
%                 ampmax = mean(max(dettrace((j-1)*16*20+1:j*16*20, 2:4)));
%                 ampmin = mean(min(dettrace((j-1)*16*20+1:j*16*20, 2:4)));
%                 midtime = dettrace((j-1)*16*20 +1 +16*20/2, 1);
%                 midtime = sprintf('%.1f',midtime);
%                 midtime = str2double(midtime);
% %                 midtime = 0.5*(dettrace((j-1)*16*20+1,1) + dettrace(j*16*20,1));
%                 objind = find(lfmaptime(:,8)==str2double(fam) & ...
%                               lfmaptime(:,9)==str2double(strcat(YEAR,JDAY)) & ...
%                               abs(lfmaptime(:,10)-midtime)< 0.5);
%                 if isempty(objind) || size(objind,1) > 1
%                     disp(midtime)
%                     disp(strcat(YEAR,JDAY))
%                     disp(objind)
%                     disp(j)
%                     break
%                 end
%                 amplf(objind,1) = ampmax;
%                 amplf(objind,2) = ampmin;
%             end
%         else
%             continue
%         end
%     end
% end
% disp('Amplitude linked to lf detections done');


%% add the above new info into the original matrix
hfnewtime = hfmaptime;
lfnewtime = lfmaptime;


%% apply distance cutoff
%%% convert to relative locations
hfrelatime = [];
lfrelatime = [];
for ifam = 1: nfam     
    tmp1 = hfnewtime(hfnewtime(:,8)==str2double(nfampool(ifam,:)),:);
    [dx, dy] = absloc2relaloc(tmp1(:,1),tmp1(:,2),loccont(ifam,1),loccont(ifam,2));
    dist = sqrt(dx.^2+dy.^2);   % distance to their own family
    hfrelatime = [hfrelatime; dx dy dist tmp1];     % now increases to 56 cols
    
    tmp2 = lfnewtime(lfnewtime(:,8)==str2double(nfampool(ifam,:)),:);
    [dx, dy] = absloc2relaloc(tmp2(:,1),tmp2(:,2),loccont(ifam,1),loccont(ifam,2));
    dist = sqrt(dx.^2+dy.^2);   % distance to their own family
    lfrelatime = [lfrelatime; dx dy dist tmp2];     % now increases to 56 cols
end

% set dist cutoff, retain detections that are within a threshold to its own family
distmaxhf = 8;
distmaxlf = 12;
% disthf = sqrt(hfrelatime(:,1).^2+hfrelatime(:,2).^2);
% distlf = sqrt(lfrelatime(:,1).^2+lfrelatime(:,2).^2);
hfrelatimecut = hfrelatime(hfrelatime(:,3)<=distmaxhf,:);
lfrelatimecut = lfrelatime(lfrelatime(:,3)<=distmaxlf,:);


% relative coordinates to fam 043
[dx, dy] = absloc2relaloc(hfrelatimecut(:,4),hfrelatimecut(:,5),-123.772167,48.493000);
hfrelatimecut = [dx dy hfrelatimecut];     % now increases to 58 cols

[dx, dy] = absloc2relaloc(lfrelatimecut(:,4),lfrelatimecut(:,5),-123.772167,48.493000);
lfrelatimecut = [dx dy lfrelatimecut];     % now increases to 58 cols

%%% save the results here
% 58 cols, format is:
%   E(002) N(002) E(own) N(own) dist(own) lon lat dep + previous 50 cols
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

fid = fopen(strcat(rstpath, '/evtloc.lzbfam.pj21.',num2str(distmaxhf),'kmdcut.',SUFFIXhf),'w+');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        hfrelatimecut');
fclose(fid);

fid = fopen(strcat(rstpath, '/evtloc.lzbfam.pj21.',num2str(distmaxlf),'kmdcut.',SUFFIXlf),'w+');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        lfrelatimecut');
fclose(fid);


%%% Also save a copy of results without dist cutoff
[dx, dy] = absloc2relaloc(hfrelatime(:,4),hfrelatime(:,5),-123.772167,48.493000);
hfrelatimenocut = [dx dy hfrelatime];

[dx, dy] = absloc2relaloc(lfrelatime(:,4),lfrelatime(:,5),-123.772167,48.493000);
lfrelatimenocut = [dx dy lfrelatime];


% format is the same as above
fid = fopen(strcat(rstpath, '/evtloc.lzbfam.pj21.nodcut.',SUFFIXhf),'w+');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        hfrelatimenocut');
fclose(fid);

fid = fopen(strcat(rstpath, '/evtloc.lzbfam.pj21.nodcut.',SUFFIXlf),'w+');
fprintf(fid,'%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %6.2f %6.2f %.4f %.4f %d %d %9.1f %10.3f %8.3f %7.2f %10.3e %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %7.3f %7.3f %7.3f %7.3f %7.3f %5.2f %5.2f %7.3f %7.3f %7.3f %6.2f %6.2f %6.2f %d %d %d %d %.4f %.4f %.4f %.4f %d %d %d %d %.3f %.3f %.3f %.3f \n',...
        lfrelatimenocut');
fclose(fid);

   








