% function synth_chaomod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THIS is to create synthetic tremor seismograms from the LFE templates 
% that have a specific distribution in time (currently is uniformly random
% in time). The composing templates have a user-specified event-size (moment)
% & frequency relation
%
%   This is the modified version of Allan's code 'chaosynth'. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2021/07/27
% Last modified date:   2021/07/27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Reads nsta(=3) template stacks. 
%Builds a synthetic seismogram for each from a user-defined window (ginput).
%shifts 2 and 3 w.r.t. 1 over a pre-defined range.
clear
clc
close all
set(0,'DefaultFigureVisible','on');
%set(0,'DefaultFigureVisible','off');

% temppath = '/Users/admin/Documents/MATLAB/';
workpath = getenv('ALLAN');
%%% choose the data path here
% datapath = workpath;
datapath = strcat(workpath,'/data-no-resp');
temppath = strcat(datapath, '/templates/PGCtrio/');


scrsz=get(0,'ScreenSize');
wid=scrsz(3)/3.1;
hite=scrsz(4);
scrat=wid/hite;

ccstack = [];
sps = 40;
templensec = 60;

fam = '002';
disp(fam);

stas=['PGC ';
      'SILB';
      'SSIB'];
nsta=size(stas,1);

for ista = 1: nsta
    fname = strcat(temppath, fam, '_', strtrim(stas(ista, :)), '_', num2str(sps), 'sps_', ...
        num2str(templensec), 's_', 'BBCCS_', 'opt_Nof_Non_Chao');
    ccstack(ista,:) = load(fname);
end
STA = ccstack;

for ista=1:nsta
    STA(ista,:)=Bandpass(STA(ista,:),sps,0.1,15,2,2,'butter');   % change 'bandpass' to 'Bandpass'
end

%%%The below aligns the templates by x-correlation
figure
hold on
plot(STA(1,:),'r')
plot(STA(2,:),'b')
plot(STA(3,:),'k')


[maxses,imaxses]=max(STA,[],2);
[minses,iminses]=min(STA,[],2);
spread=maxses-minses;
% zcrosses=round(0.5*(imaxses+iminses));  % rough, assuming symmetry, Chao 2021/07/16
zcrosses = zeros(nsta,1);
for ista = 1:nsta
    seg = STA(ista, iminses(ista): imaxses(ista));  % for zero-crossing timing, only use the main station
    [~,zcrosses(ista)] = min(abs(seg));
    zcrosses(ista) = zcrosses(ista)-1+iminses(ista);  % convert to global index
end
sampbef=6*sps;
sampaft=10*sps;
is=zcrosses-sampbef;
ie=zcrosses+sampaft;
for ista=1:nsta
    STAtmp(ista,:)=STA(ista,is(ista):ie(ista));
end
mshift=10;
tempxc(1,:)=xcorr(STAtmp(2,:),STAtmp(3,:),mshift,'coeff');
tempxc(2,:)=xcorr(STAtmp(1,:),STAtmp(2,:),mshift,'coeff'); %shift STAtmp(3,:) to right for positive values
tempxc(3,:)=xcorr(STAtmp(1,:),STAtmp(3,:),mshift,'coeff'); %shift STAtmp(2,:) to right for positive values
[~,imax]=max(tempxc,[],2);
imax=imax-(mshift+1) %This would produce a slightly different shift, if filtered seisms were used.
for ista=2:nsta
    STAtmp(ista,mshift+1:end-(mshift+1))=STAtmp(ista,mshift+1-imax(ista):end-(mshift+1)-imax(ista));
end
for ista=1:nsta
    STAtmp(ista,:)=STAtmp(ista,:)/spread(ista);
end
figure
hold on
plot(STAtmp(1,:),'r')
plot(STAtmp(2,:),'b')
plot(STAtmp(3,:),'k')
%%%The above aligns the templates by x-correlation

%%%The below determines the template window, and makes sure it has zero mean with tapered ends
% for tapering, could also use 'tukeywin', Chao 2021/07/16
% greenlen is its length
%according analysis to the template, the zero-crossing to the end of coda is safe to set as 8 s
zcross = 241;
greenlen = pow2(9);
tlim = [zcross+8*sps-greenlen+1  zcross+8*sps];
tgreen = zcross-tlim(1)+1;  % zero-crossing time of the final truncated tempalte
% tlim=[180 899]; %just for consistency; could go back to picking
% greenlen=tlim(2)-tlim(1)+1;
green=zeros(nsta,greenlen);
smooth=sps/2;
for ista=1:nsta
    green(ista,:)=STAtmp(ista,tlim(1):tlim(2));
    green(ista,:)=green(ista,:)-mean(green(ista,:));    % remove mean
    x=(0:pi/(2*smooth):pi/2-pi/(2*smooth));    
    green(ista,1:smooth)=sin(x).*green(ista,1:smooth);  % left taper
    x=fliplr(x);
    green(ista,end-(smooth-1):end)=sin(x).*green(ista,end-(smooth-1):end);  % right taper
    green(ista,:)=green(ista,:)-mean(green(ista,:));    % remove mean
    green(ista,:)=green(ista,:)/max(abs(green(ista,:)));    % normalize
end
mean(green,2)
subplot(2,1,1)
hold on
plot(green(1,:),'r')
plot(green(2,:),'b')
plot(green(3,:),'k')
mx=max([abs(green(1,:)) abs(green(2,:))]);
xlim([0 greenlen])
ylim([-mx mx])
box on

for ista=1:nsta
    greenf(ista,:)=Bandpass(green(ista,:),100,1.8,18,2,2,'butter');   % change 'bandpass' to 'Bandpass'
    greenf(ista,:)=greenf(ista,:)-mean(greenf(ista,:)); % remove mean again for caution
end

subplot(2,1,2)
hold on
plot(greenf(1,:),'r')
plot(greenf(2,:),'b')
plot(greenf(3,:),'k')
mx=max([abs(greenf(1,:)) abs(greenf(2,:))]);
%%%The above determines the template window, and makes sure it has zero mean with tapered ends

%%%The below cross-correlates the templates (Green's functions), just FYI
% running CC using a window length of 'cclen', Chao 2021/07/26
    cclen=sps/2;
    samples=1:length(greenf(1,:));
    ST1auto=greenf(1,:).*greenf(1,:);
    ST1sq=cumsum(ST1auto);
    ST2auto=greenf(2,:).*greenf(2,:);
    ST2sq=cumsum(ST2auto);
    ST3auto=greenf(3,:).*greenf(3,:);
    ST3sq=cumsum(ST3auto);
    STA12tr=greenf(1,:).*greenf(2,:);            
    STA13tr=greenf(1,:).*greenf(3,:);
    STA23tr=greenf(2,:).*greenf(3,:);
    cumsumSTA12tr=cumsum(STA12tr);
    cumsumSTA13tr=cumsum(STA13tr);
    cumsumSTA23tr=cumsum(STA23tr);
    ST12num=cumsumSTA12tr(cclen+1:length(greenf(1,:)))-cumsumSTA12tr(1:length(greenf(1,:))-cclen);
    ST13num=cumsumSTA13tr(cclen+1:length(greenf(1,:)))-cumsumSTA13tr(1:length(greenf(1,:))-cclen);
    ST23num=cumsumSTA23tr(cclen+1:length(greenf(1,:)))-cumsumSTA23tr(1:length(greenf(1,:))-cclen);
    ST1den=ST1sq(cclen+1:length(greenf(1,:)))-ST1sq(1:length(greenf(1,:))-cclen);
    ST2den=ST2sq(cclen+1:length(greenf(1,:)))-ST2sq(1:length(greenf(1,:))-cclen);
    ST3den=ST3sq(cclen+1:length(greenf(1,:)))-ST3sq(1:length(greenf(1,:))-cclen);
    ST12n=ST12num./realsqrt(ST1den.*ST2den);
    ST13n=ST13num./realsqrt(ST1den.*ST3den);
    ST23n=ST23num./realsqrt(ST2den.*ST3den);
    alln=(ST12n+ST13n+ST23n)/3;
    %alln(alln<0)=-10^4*yma; %just so they don't plot.
    plot(samples(cclen/2+1:greenlen-cclen/2),mx*alln,'co','markersize',2); % scale with data amp.
xlim([0 greenlen])
ylim([-mx mx])
box on
xc23=xcorr(greenf(2,:),greenf(3,:),10,'coeff');
xc13=xcorr(greenf(1,:),greenf(3,:),10,'coeff');
xc12=xcorr(greenf(1,:),greenf(2,:),10,'coeff');
[ccmax23,imax23]=max(xc23)
[ccmax13,imax13]=max(xc13)
[ccmax12,imax12]=max(xc12)
%%%The above cross-correlates the template (Green's functions), just FYI

%%
%%%Choose the distribution of event magnitudes. I'm not sending you most of these.  Power-law,
% log-normal, uniform (every event the same size), "discrete Ide", exponential.
distr='PL'
%distr='LN'
% distr='UN'
%distr='DI'
%distr='EX'
%print('-depsc',['b=',int2str(b),'_template.eps'])
% not sure about '+3', 'greenlen/sps' is for deleting, --Chao
Twin=0.5*3600+3+2*ceil(greenlen/sps); %Early and late portions will be deleted. Twin includes only the early portion. In seconds.
winlen=Twin*sps+1;  
skiplen=greenlen;   % length that will be skipped in output, --Chao
%How saturated the seismogram:  imagine the full win is covered by non-overlapping templates --Chao
satn=2.5*Twin %a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
    %Twin is window duration in seconds. Events can fall within Twin of
    %the start, but the synthetics will go to Twin*(sample rate)+Greenlen to
    %avoid checking for subscript overrrun.  When "synths" are written from "synth" in the
    %subroutine, Greenlen from the start and end will not be written.
mjig=0; %5; %max number of samples to jigger 2 and 3 w.r.t. 1.  AT 100 SPS!  Currently not used
fracelsew=0.; %0.25 %0.5; %0.6; %The fraction of elsewhere events to local events.  "Elsewhere events" have random times at EACH station.
%writes=[1e7]; %Just a single entry for Discrete Ide.  Value not relevant.
nsat=[0.1 0.4 1 2 4 10 100]; %One catalog at each of these levels of saturation.
writes=round(nsat*satn) %in degree of saturation, i.e. the number of composing templates --chao
seed=randi(100,1);
rng(seed,'twister');
% rng('default');
% with some extra length for deleting and shifting of 2nd and 3rd sta
% synth=2.e-7*(rand(nsta,winlen+greenlen+2*10)-0.5); %a bit of noise for later cross-correlation of under-saturated seismograms.
synth=2.e-3*(rand(nsta,winlen+greenlen+2*10)-0.5); %a bit of noise for later cross-correlation of under-saturated seismograms.
%+2*10 a little extra, for jiggering 2nd & 3rd stations.% for ista=1:nsta

nouts=length(writes); %how many output files
% seed=round(writes(1)/5e3); %for random number generator
if strcmp(distr,'PL') || strcmp(distr,'UN') %b>150 for uniform %Go here for a power-law size distribution or uniform sizes
    b=999. %>150 for uniform size distribution
%    xygrid=load('xygridArmb'); %Made from /ARMMAP/MAPS/2021/interpgrid.m; format PGSS, PGSI, dx, dy
%    size(xygrid)
%    %sourcestyle='ellipse';
%    sourcestyle='rectangle';
%    if strcmp(sourcestyle,'rectangle') %The following rotates 45˚ and looks for -5 < x' < 1 and 2.25 < y' < 3.5
%        angrot=-45*pi/180;
%        loc=complex(xygrid(:,3),xygrid(:,4));
%        locrot=loc*exp(1i*angrot);
%        xygrid(:,5)=real(locrot);
%        xygrid(:,6)=imag(locrot);
%        ll=[-5 2.25]; ur=[1 3.5];
%        xygrid(xygrid(:,5)<ll(1) | xygrid(:,5)>ur(1) | xygrid(:,6)<ll(2) | xygrid(:,6)>ur(2),:)=[]; %This limits xygrid to the desired range.  PGSS, PGSI, x, y, x', y'
%    elseif strcmp(sourcestyle,'ellipse') %The following adds 0.5 km to y, rotates 45˚, and looks for an elliptical region with major semi-axis y' < 3 km and minor semi-axis < 2 km
%        angrot=-45*pi/180;
%        shiftor=[0 -0.5]; %(in km)
%        loc=complex(xygrid(:,3)-shiftor(1),xygrid(:,4)-shiftor(2));
%        locrot=loc*exp(1i*angrot);
%        xygrid(:,5)=real(locrot);
%        xygrid(:,6)=imag(locrot);
%        xaxis=0.5; %(semi-, in km)
%        yaxis=0.5; %(semi-, in km)
%        xygrid(sqrt((xygrid(:,5)/xaxis).^2 + (xygrid(:,6)/yaxis).^2) > 1,:)=[]; %This limits xygrid to the desired range.  PGSS, PGSI, x, y, x', y'
%    end
%    fid = fopen('tmp.grd','w');
%    fprintf(fid,'%8.3f %8.3f %7.2f %7.2f %8.3f %8.3f\n',xygrid');
%    fclose(fid);
%    size(xygrid)
%     [synths,mommax]=plaw3b(writes,winlen,skiplen,synth,green,b,mjig,fracelsew,seed); 
    [synths,mommax,greensts]=plaw3b(writes,winlen,skiplen,synth,greenf,b,mjig,fracelsew,seed); 
end
if strcmp(distr,'LN')
    sig=1.;
    [synths,mommax]=lognorm(writes,winlen,skiplen,synth,green,sig,mjig,fracelsew,seed); 
end
if strcmp(distr,'DI')
    stepsize=5; %Sets the step size for the random walk
    frac=1; %Sets the s.d. of the variation as frac*sqrt(moment rate)
    desiredmax=3000; %2000; %Used to shift the center of the Gaussian, to limit large moment rates. 
    alpha=1; %Center of the Gaussian is shifted as -sd*(STF(i-1)/desiredmax)^alpha
    av=30*sps; %in samples, the length of time over which the average moment rate is determined
    
    stepsize=5; %Sets the step size for the random walk
    frac=1; %Sets the s.d. of the variation as frac*sqrt(moment rate)
    desiredmax=1e9; %2000; %Used to shift the center of the Gaussian, to limit large moment rates. 
    alpha=1; %Center of the Gaussian is shifted as -sd*(STF(i-1)/desiredmax)^alpha
    av=1; %in samples, the length of time over which the average moment rate is determined
    
    [synths,cumSTF]=DIde(writes,winlen,skiplen,synth,stepsize,desiredmax,frac,alpha,av,green,mjig,seed);
end
if strcmp(distr,'DI')
    for inwrites=1:length(writes)
        n=cumSTF;
        fid = fopen(['3STAS/STAS.',distr,'.dt',int2str(stepsize),'.dmax',int2str(desiredmax),'.av',int2str(av),'.alp',num2str(alpha),'.frac',num2str(frac),'.',int2str(sps),'sps.mjig',int2str(mjig),'else',num2str(fracelsew,2),'_N10e',num2str(log10(n),3)],'w');
        towrite=squeeze(synths(:,inwrites,:));
        fprintf(fid,'%13.5e %13.5e %13.5e\n', towrite);
        fclose(fid);
    end
end
if strcmp(distr,'LN')
    for inwrites=1:length(writes)
        n=writes(inwrites);
        fid = fopen(['3STAS/STAS.',distr,'.sig',num2str(sig,2),'.',int2str(sps),'sps.mjig',int2str(mjig),'else',num2str(fracelsew,2),'_N10e',num2str(log10(n),3)],'w');
        towrite=squeeze(synths(:,inwrites,:));
        fprintf(fid,'%13.5e %13.5e %13.5e\n', towrite);
        fclose(fid);
    end
end
if strcmp(distr,'UN')
    for inwrites=1:length(writes)
        n=writes(inwrites);
        if strcmp(sourcestyle,'ellipse')
            fid = fopen(['3STAS/STAS.',distr,'.',int2str(sps),'sps.',sourcestyle(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites))],'w');
        elseif strcmp(sourcestyle,'rectangle')
            fid = fopen(['3STAS/STAS.',distr,'.',int2str(sps),'sps.',sourcestyle(1:3),'_',num2str(ll(1)),num2str(ll(2)),num2str(ur(1)),num2str(ur(2)),'.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites))],'w');
        end
        towrite=squeeze(synths(:,inwrites,:));
        fprintf(fid,'%13.5e %13.5e %13.5e\n', towrite);
        fclose(fid);
        if strcmp(sourcestyle,'ellipse')
            fid = fopen(['3STAS/STAS.',distr,'.',int2str(sps),'sps.',sourcestyle(1:3),'_',num2str(xaxis),'-',num2str(yaxis),'.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'_sources'],'w');
        elseif strcmp(sourcestyle,'rectangle')
            fid = fopen(['3STAS/STAS.',distr,'.',int2str(sps),'sps.',sourcestyle(1:3),'_',num2str(ll(1)),num2str(ll,2),num2str(ur(1)),num2str(ur,2),'.else',num2str(fracelsew,2),'nsat',num2str(nsat(inwrites)),'_sources'],'w');
        end
        truncsources=sortrows(sources(1:n,:),1);
        fprintf(fid,'%13i %5i %4i\n', truncsources');
        fclose(fid);
    end
end

%%
figure
nrow = length(writes)+1;
ncol = 1;
subplot(nrow,ncol,1)
hold on
tmp = synth;
plot(tmp(1,:),'r')
plot(tmp(2,:),'b')
plot(tmp(3,:),'k')
xlim([0 0.4e4])

for i = 1: length(writes)
subplot(nrow,ncol,i+1)
hold on
tmp = synths(:,i,:);
plot(tmp(1,:),'r')
plot(tmp(2,:),'b')
plot(tmp(3,:),'k')
xlim([0 0.2e4])
end

%% 
for i = 1: length(writes)
    tmp = squeeze(synths(:,i,:));
    ampmax = median(max(abs(tmp),[],2));
    ntempsqrt = sqrt(writes(i));
    ntempsqrt/ampmax
end    
    
%% extract and validate the added impulses of template 
insat = 1;
wlen = 20;
ista = 2; 

indstsig = 1;  % starting index of the simulated signal to test  
indedsig = wlen*sps+indstsig-1; % ending index of the simulated signal to test  
greenst = sort(greensts{insat})'; % the starting index of each added template, context of full length
%you don't need all impulses, only some of them contribute to the length of truncated record 
greenst = greenst(greenst>=skiplen-greenlen+indstsig & greenst<=indedsig+skiplen-1);   
greenzc = greenst+tgreen; % index of approximate zero-crossing
impamp = zeros(max(greenzc),1);
for i = 1: length(greenzc)
  impamp(greenzc(i)) = impamp(greenzc(i))+1;
end

%note that some length of the full simulation is skipped
sigplnoi = squeeze(synths(:,insat,:));  % simulated signal with white noise
sig = sigplnoi-synth(:,skiplen:winlen); % subtract the white noise

sig = sig(ista,indstsig:indedsig)';
sigconv = conv(greenf(ista,:)',impamp,'full');
indtrc = tgreen+skiplen;  % starting index for truncation
sigconv = sigconv(indtrc+indstsig-1:indedsig+indtrc-1);  % cut accordingly

% figure
% plot(sig,'k'); hold on
% plot(sigconv+3,'b');
% plot(sigconv-sig+1,'r-')

figure
subplot(311)
plot(greenf(ista,:),'r');
xlim([0 greenlen]);
legend('Template');

subplot(312)
imptemp = find(impamp>0);
p1=stem(imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
% p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
ax = gca;
plot([indstsig indstsig],ax.YLim,'--','color',[.5 .5 .5]);
plot([indedsig indedsig],ax.YLim,'--','color',[.5 .5 .5]);
legend(p1,'Synthetic random impulses');

subplot(313);
sigplnoi = sigplnoi(1,indstsig:indedsig)';
plot(sig,'b'); hold on
plot(sigconv,'k');
legend('Truncated synthetic signal','Truncated signal from convolution');

%% THIS section is only for grouping meeting presentation
%%%For showing purpose
widin = 6;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
nrow = 4;
ncol = 1;
f = initfig(widin,htin,nrow,ncol);  % create a figure with subplots

pltxran = [0.1 0.9]; pltyran = [0.1 0.9];
pltxsep = 0.02; pltysep = 0.04; 
%get the locations for each axis
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax = f.ax(1);
plot(ax,greenf(ista,:),'r'); hold on
% scatter(tgreen,greenf(ista,tgreen),20,'k','filled');
% xlim([0 greenlen]);
xlim(ax,[0 greenlen]);
nolabels(ax,1);
legend(ax,'Template');
ylabel(ax,'Amplitude');
shrink(ax,wlen*sps/greenlen,1);
ax.Position(1)=axpos(1,1);
longticks(ax,1);
axsym(ax,2);
axranexp(ax,6,10);

ax = f.ax(2);
imptemp = find(impamp>0);
p1=stem(ax,imptemp-skiplen, impamp(imptemp),'b','MarkerSize',4); hold on;
% p1=stem((1:size(impamp))-skiplen, impamp,'b','MarkerSize',4); hold on;
% ax = gca;
% plot([indstsig indstsig],ax.YLim,'--','color',[.5 .5 .5]);
% plot([indedsig indedsig],ax.YLim,'--','color',[.5 .5 .5]);
nolabels(ax,1);
longticks(ax,2);
legend(ax,p1,'Synthetic random impulses');
xlim(ax,[0 wlen*sps]);
ylabel(ax,'Amplitude');
axranexp(ax,2,10);

ax = f.ax(3);
% plot(sig,'b'); hold on
plot(ax,sigconv,'k');
% legend('Truncated synthetic signal','Truncated signal from convolution');
nolabels(ax,1);
longticks(ax,2);
legend(ax,'Truncated signal from convolution');
ylabel(ax,'Amplitude');
axranexp(ax,6,10);

noistd = 5.e-2;
rng(seed,'twister');
noi = noistd*randn(length(sigconv),1);
sigplnoi = sigconv+noi;
ax = f.ax(4);
% plot(sig,'b'); hold on
plot(ax,sigplnoi,'k');
% legend('Truncated synthetic signal','Truncated signal from convolution');
longticks(ax,2);
legend(ax,sprintf('Gaussian noise with std=%.2f added',noistd));
xlabel(ax,sprintf('Samples at %d Hz',sps));
ylabel(ax,'Amplitude');
axranexp(ax,6,10);
axranexp(ax,6,10);



wlet = greenf(ista,:)';
sig = sigplnoi;
% noi = [];

%detrend and taper
sig = detrend(sig);
fractap = 0.05; % if fractap is >=1, n-point von Hann window is returned
ptstap = fractap/2*size(sig,1); % if fractap is >=1, n-point von Hann window is returned
w = tukeywin(size(sig,1),fractap);
% sig = w.* sig;
%detrend again for caution
sig=detrend(sig);


dt = 1/sps;  % sampling interval
twlet = tgreen*dt;
width = 2.5;  % width for Gaussian filter
dres_min = 0.5;  % tolerance, percentage change in residual per iteration
mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
nit_max = 20;  % max numer of iterations 
npul_max = round(2.5*wlen*nsat(insat));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
fpltit = 0;  % plot flag for each iteration
fpltend = 1;  % plot flag for the final iteration
fcheck = 0;
rcc = [];  % running CC between diff. stations
rcc = ones(length(sig) ,1); % for testing purpose    
[sigdecon,pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  iterdecon(sig,wlet,rcc,noi,dt,twlet,width,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fcheck);

%%
widin = 6;  % maximum width allowed is 8.5 inches
htin = 5;   % maximum height allowed is 11 inches
nrow = 4;
ncol = 1;
f = initfig(widin,htin,nrow,ncol);  % create a figure with subplots

pltxran = [0.1 0.9]; pltyran = [0.1 0.9];
pltxsep = 0.02; pltysep = 0.04; 
%get the locations for each axis
axpos = optaxpos(f,nrow,ncol,pltxran,pltyran,pltxsep,pltysep);

ax = f.ax(1);
plot(ax,greenf(ista,:),'r'); hold(ax, 'on');
% scatter(tgreen,greenf(ista,tgreen),20,'k','filled');
% xlim([0 greenlen]);
xlim(ax,[0 greenlen]);
nolabels(ax,1);
legend(ax,'Template');
ylabel(ax,'Amplitude');
shrink(ax,wlen*sps/greenlen,1);
ax.Position(1)=axpos(1,1);
longticks(ax,1);
axsym(ax,2);
axranexp(ax,6,10);
nolabels(ax,1);
hold(ax, 'off');

ax = f.ax(2);
imptemp = find(impamp>0);
stem(ax,imptemp-skiplen, impamp(imptemp),'k','MarkerSize',4); hold(ax, 'on');
imptemp = find(sigdecon>0); 
stem(ax,imptemp,sigdecon(imptemp),'b','MarkerSize',4);
legend(ax,'Synthetic random impulses','Deconvolved impulses');
xlim(ax,[0 wlen*sps]);
ylabel(ax,'Amplitude');
longticks(ax,2);
axranexp(ax,2,10);
nolabels(ax,1);
hold(ax, 'off');

ax = f.ax(3);
plot(ax,sig,'k'); hold(ax, 'on');
plot(ax,pred, 'b');
xlim(ax,[0 wlen*sps]);
legend(ax,'Signal','Prediction');
ylabel(ax,'Amplitude');
longticks(ax,2);
axsym(ax,2);
axranexp(ax,6,10);
nolabels(ax,1);
hold(ax, 'off');

ax = f.ax(4);
p1=plot(ax,res, 'c-'); ylim(ax.YLim); hold(ax, 'on');
p2=plot(ax,noi,'-','color',[.7 .7 .7]);
xlim(ax,[0 wlen*sps]);
legend(ax,[p1,p2],'Residual','Noise');
xlabel(ax,sprintf('Samples at %d Hz',sps));
ylabel(ax,'Amplitude');
longticks(ax,2);
ylim(ax,f.ax(3).YLim);
hold(ax, 'off');

%% 
% sig = sig;
% noi = [];
% sig = sigplnoi;
% noi = synth(1,skiplen:winlen)';
%detrend and taper
sig = detrend(sig);
fractap = 0.05; % if fractap is >=1, n-point von Hann window is returned
ptstap = fractap/2*size(sig,1); % if fractap is >=1, n-point von Hann window is returned
w = tukeywin(size(sig,1),fractap);
sig = w.* sig;
%detrend again for caution
sig=detrend(sig);
  
wlet = greenf(ista,:)';

% sig = zeros(size(sig)); % for testing purpose
% sig = sig+0.6*wlet(1:size(sig));
% sig = sig+0.8*wlet(41:size(sig)+40);
% sig = sig+1*wlet(101:size(sig)+100);

% sig(1:1+length(wlet)-100-1) = sig(1:1+length(wlet)-100-1)+0.6*wlet(101:end);
% sig(800:800+length(wlet)-1) = sig(800:800+length(wlet)-1)+0.8*wlet;
% sig(400:400+length(wlet)-1) = sig(400:400+length(wlet)-1)+0.4*wlet;
% sig(1601:end) = sig(1601:end)+1*wlet(1:400);


%% try to deconvolve it, using iterative time-domain deconvolution method
dt = 1/sps;  % sampling interval
twlet = tgreen*dt;
width = 2.5;  % width for Gaussian filter
dres_min = 0.5;  % tolerance, percentage change in residual per iteration
mfit_min = 5e-1;  % tolerance, norm of the residual/misfit
nit_max = 10000;  % max numer of iterations 
nimp_max = round(2.5*wlen*nsat(insat));%a single peak is ~20 samples wide; maybe a little less (at 100 sps). ~0.4s, 1/0.4=2.5
fpltit = 1;  % plot flag for each iteration
fpltend = 1;  % plot flag for the final iteration
rcc = [];  % running CC between diff. stations
rcc = ones(length(sig) ,1); % for testing purpose    
[sigdecon,pred,res,dresit,mfitit,ampit,nit,fighdl] = ...
  iterdecon(sig,wlet,rcc,noi,dt,twlet,width,dres_min,mfit_min,nit_max,npul_max,fpltit,fpltend,fcheck);

[coef,lag] = xcorr(impamp(skiplen+indstsig:skiplen+indedsig), sigdecon, sps/4, 'coeff');
[mcoef, ind] = max(coef);
lagsamp = lag(ind);
mcoef
lagsamp

ax = fighdl{2}.ax(2);
hold(ax,'on');
p1=stem(ax,impamp(skiplen+indstsig:skiplen+indedsig), 'k','MarkerSize',4); 
text(ax,0.05,0.9,sprintf('lag: %d, max coef: %.2f', lagsamp,mcoef),'Units','normalized',...
  'HorizontalAlignment','left');
legend(ax,'Modeled','True');
hold(ax,'off');


%% try the various spectral division methods, one-time division in freq-domain
dt = 1/sps;  % sampling interval
twlet = tgreen*dt;
stdec1 = specdiv_damp(sig,wlet,noi,dt,twlet);
figure
subplot(311)
stem(impamp(skiplen+indstsig:skiplen+indedsig),'k'); hold on
stem(sigdecon,'b');
stem(stdec1,'r');
legend('True','Iterdecon','Specdiv\_damp');

stdec2 = specdiv_water(sig,wlet,noi,dt,twlet);
subplot(312)
stem(impamp(skiplen+indstsig:skiplen+indedsig),'k'); hold on
stem(sigdecon,'b'); hold on
p1=stem(stdec2,'r');
legend(p1,{'Specdiv\_water'});

stdec3 = specdiv_arraycond(sig,wlet,noi,dt,twlet);
subplot(313)
stem(impamp(skiplen+indstsig:skiplen+indedsig),'k'); hold on
stem(sigdecon,'b'); hold on
p1=stem(stdec3,'r');
legend(p1,'Specdiv\_arraycond');









