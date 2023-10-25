function [STAEuse,STANuse,STAUuse,STAtuse,STAEpsig,STANpsig,STAUpsig,STAtpsig,...
          STAEssig,STANssig,STAUssig,STAtssig,STAEnoi,STANnoi,STAUnoi,STAtnoi] = ...
    rdhinetevt(sacpath,stainfo,evt,sps,lo,hi,npo,npa,wlensec,wlensecbf,pltflag,velmod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to read the regular earthaquake event data from Japan Hinet
% (after or before instrument response removal), bandpass filtering, resample if
% needed. Return the 120-s long window of E-N-U components, p/S signal window, 
% s/S signal window, whose length is as requested
%
% 
%   NOTICE:
%       1. The time zone now used in both tremor, regular events catalog and the Hinet
%           data are JST (UT+9), Japan Standard Time = UTC + 9 h, which are consistent 
%           to each other. Unless specified, conversion to UTC frame is unnecessary
%
%
%
%
% By Chao Song, chaosong@princeton.edu
% First created date:   2020/03/05
% Last modified date:   2020/03/05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% pre-define the output, incase the data is bad somehow
%%% this station can't be used, moreover, if this station isessential, then this event should be
%%% skipped
tlensec = 120;
tlen = tlensec*sps;

% STAEuse = nan(tlen,1);
% STANuse = nan(tlen,1);
% STAUuse = nan(tlen,1);
% STAtuse = nan(tlen,1);
% STAEnoi = nan(wlensec*sps,1);
% STANnoi = nan(wlensec*sps,1);
% STAUnoi = nan(wlensec*sps,1);
% STAtnoi = nan(wlensec*sps,1);
% STAEpsig = nan(wlensec*sps,1);
% STANpsig = nan(wlensec*sps,1);
% STAUpsig = nan(wlensec*sps,1);
% STAtpsig = nan(wlensec*sps,1);
% STAEssig = nan(wlensec*sps,1);
% STANssig = nan(wlensec*sps,1);
% STAUssig = nan(wlensec*sps,1);
% STAtssig = nan(wlensec*sps,1);
STAEuse = [];
STANuse = [];
STAUuse = [];
STAtuse = [];
STAEnoi = [];
STANnoi = [];
STAUnoi = [];
STAtnoi = [];
STAEpsig = [];
STANpsig = [];
STAUpsig = [];
STAtpsig = [];
STAEssig = [];
STANssig = [];
STAUssig = [];
STAtssig = [];

% %%% DEFINE a quality flag, if flag is 1, data is good, if 0, data is bad in some way
% 
% qualflag = 1;
        
stnm = char(stainfo{1,1});    % station name, string
stlo = stainfo{1,2};    % station longitude
stla = stainfo{1,3};    % station latitude
    
for i = 1: size(stnm,1)
% i=1;

    clear tmpe tmpn tmpu tmpef tmpnf tmpuf tmpefd tmpnfd tmpufd
    datestr = num2str(evt(1));
    dirnm = datestr(3:end);
%     staefnm = strcat(stnm(i,:),'.E');
%     stanfnm = strcat(stnm(i,:),'.N');
%     staufnm = strcat(stnm(i,:),'.U');
    staefnm = strcat(stnm(i,:),'.E.new');
    stanfnm = strcat(stnm(i,:),'.N.new');
    staufnm = strcat(stnm(i,:),'.U.new');
    staepath = fullfile(sacpath,dirnm,staefnm);
    stanpath = fullfile(sacpath,dirnm,stanfnm);
    staupath = fullfile(sacpath,dirnm,staufnm);
    
    if isfile(staepath)
        [tmpe,hdr,~,~,time]=readsac(staepath,0,'l');   %%% Function 'readsac'
        [tmpn,~,~,~,~]=readsac(stanpath,0,'l');
        [tmpu,~,~,~,~]=readsac(staupath,0,'l');
    else
        return
    end
    
%     sps = 50;
    spstrue = round(1./hdr.DELTA);

    %%% ALWAYS filter first, then resample, to avoid time-domain waveform distortion
    % filtering
%     lo = 0.05;
%     hi = 20;
%     npo = 2;
%     npa = 2;
    tmpef = Bandpass(tmpe,spstrue,lo,hi,npo,npa,'butter');
    tmpnf = Bandpass(tmpn,spstrue,lo,hi,npo,npa,'butter');
    tmpuf = Bandpass(tmpu,spstrue,lo,hi,npo,npa,'butter');
   
    % resample if needed
    [num, denom] = rat(sps/spstrue);
    tmpefd = resample(tmpef,num,denom);   % times = num/denom
    tmpnfd = resample(tmpnf,num,denom);
    tmpufd = resample(tmpuf,num,denom);
    time = linspace(round(time(1)), round(time(end)),...
                           round(length(time)*sps/spstrue))';
                           
    STAE(:,i) = tmpefd;     % 1-day E component data
    STAN(:,i) = tmpnfd;     % 1-day N component data
    STAU(:,i) = tmpufd;     % 1-day U component data
    STAt(:,i) = time;       % 1-day time data
            
%     tmped = resample(tmpe,num,denom);
%     tmpedf = Bandpass(tmped,sps,lo,hi,npo,npa,'butter');  
%     figure
%     hold on
%     plot(tmpefd,'r');
%         plot(tmpedf,'b');
%     figure
%     plot(tmpefd-tmpedf,'k');

end

dt = 1/sps;
otsec = evt(1,5)*3600+evt(1,6)*60+evt(1,7);     % origin time
otind = find(abs(otsec-time)<=dt,1,'first');    % origin time index
if isempty(otind)
    return
end

STAEuse = STAE(otind: min(otind+tlen-1, size(time,1)), :);
STANuse = STAN(otind: min(otind+tlen-1, size(time,1)), :);
STAUuse = STAU(otind: min(otind+tlen-1, size(time,1)), :);
STAtuse = STAt(otind: min(otind+tlen-1, size(time,1)), :);
if size(STAEuse,1) < tlen
    return
end

wlensecaf = wlensec-wlensecbf;
% tmpstind = max(otind-wlensecaf*sps,1);
STAEnoi = STAE(otind-wlensecaf*sps: min(otind+wlensecbf*sps-1, size(time,1)), :);
STANnoi = STAN(otind-wlensecaf*sps: min(otind+wlensecbf*sps-1, size(time,1)), :);
STAUnoi = STAU(otind-wlensecaf*sps: min(otind+wlensecbf*sps-1, size(time,1)), :);
STAtnoi = STAt(otind-wlensecaf*sps: min(otind+wlensecbf*sps-1, size(time,1)), :);
if size(STAEuse,1) < tlen
    return
end

% several available velocity models
if isequal(velmod,'jma2001')
    vmod = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/jma2001_sparse.taup';
elseif isequal(velmod,'ukawa83')
    vmod = '/home/data2/chaosong/Seisbasics/TauP-2.4.5/StdModels/ukawa83.taup';
elseif isequal(velmod,'ak135')
    vmod = 'ak135';
end

for i = 1: size(stnm,1)
    % calculate travel time                      
    tt=tauptime('mod',vmod,'dep',evt(1,10),'ph','p','evt',[evt(1,9) evt(1,8)],...
                'sta',[stla(i) stlo(i)]);
    if isempty(tt)
        fprintf('No p at station %s, using P instead. \n',stnm(i,:));
        tt=tauptime('mod',vmod,'dep',evt(1,10),'ph','P','evt',[evt(1,9) evt(1,8)],...
                    'sta',[stla(i) stlo(i)]);
    end
    ttp = tt(1).time;   % travel time p or P

    tt=tauptime('mod',vmod,'dep',evt(1,10),'ph','s','evt',[evt(1,9) evt(1,8)],...
                'sta',[stla(i) stlo(i)]);
    if isempty(tt)
        fprintf('No s at station %s, using S instead. \n',stnm(i,:));
        tt=tauptime('mod',vmod,'dep',evt(1,10),'ph','S','evt',[evt(1,9) evt(1,8)],...
                    'sta',[stla(i) stlo(i)]);
    end
    tts = tt.time;   % travel time s
    
    atpsec = otsec+ttp;    % arrival time p
    atssec = otsec+tts;    % arrival time s
    
    atpind(i) = find(abs(atpsec-time)<=dt,1,'first');    % origin time index
    atsind(i) = find(abs(atssec-time)<=dt,1,'first');
    
    if isempty(atpind(i)) || isempty(atsind(i))
        disp('Bad event! P or S arrival exceeds the event day')
        return
    end
    
    % window contains EQ p arrival + tremor + noise
    if atpind(i)+wlensecaf*sps-1 <= size(time,1)
        STAEpsig(:,i) = STAE(atpind(i)-wlensecbf*sps: atpind(i)+wlensecaf*sps-1, i);
        STANpsig(:,i) = STAN(atpind(i)-wlensecbf*sps: atpind(i)+wlensecaf*sps-1, i);
        STAUpsig(:,i) = STAU(atpind(i)-wlensecbf*sps: atpind(i)+wlensecaf*sps-1, i);
        STAtpsig(:,i) = STAt(atpind(i)-wlensecbf*sps: atpind(i)+wlensecaf*sps-1, i);
    else
        disp('Bad event! P arrival signal is shorter than needed')
        return
    end
    
    % window contains EQ s arrival + tremor + noise
    if atsind(i)+wlensecaf*sps-1 <= size(time,1)
        STAEssig(:,i) = STAE(atsind(i)-wlensecbf*sps: atsind(i)+wlensecaf*sps-1, i);
        STANssig(:,i) = STAN(atsind(i)-wlensecbf*sps: atsind(i)+wlensecaf*sps-1, i);
        STAUssig(:,i) = STAU(atsind(i)-wlensecbf*sps: atsind(i)+wlensecaf*sps-1, i);
        STAtssig(:,i) = STAt(atsind(i)-wlensecbf*sps: atsind(i)+wlensecaf*sps-1, i);
    else
        disp('Bad event! S arrival signal is shorter than needed')
        return
    end
   
end


if pltflag && ~isempty(STAEuse) && ~isempty(atpind) && ~isempty(atsind)
    
    [scrsz, res] = pixelperinch(1);

    figure
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 11;   % maximum height allowed is 11 inches
    set(gcf,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    
    subplot(311)
    hold on
    for i = 1: size(stnm,1)
        % plot waveform
        plot(STAEuse(:,i)/max(STAEuse(:,i))+i,'color',[0.5 0.5 0.5]);
        % plot p wave arrival
        plot([atpind(i)-otind atpind(i)-otind],[-1+i 1+i],'r--','linew',1);
        text(atpind(i)-otind - 300, i+0.5, 'p/P','fontsize',10,'color','r');
        % indicate the p wave signal window
        plot([atpind(i)-otind-wlensecbf*sps, atpind(i)-otind+wlensecaf*sps-1],[i+0.5 i+0.5], ...
             'k-','linew',3);
        % plot s wave arrival
        plot([atsind(i)-otind atsind(i)-otind],[-1+i 1+i],'b--','linew',1);
        text(atsind(i)-otind - 300, i+0.5, 's/S','fontsize',10,'color','b');
        % indicate the s wave signal window
        plot([atsind(i)-otind-wlensecbf*sps, atsind(i)-otind+wlensecaf*sps-1],[i+0.5 i+0.5], ...
             'k-','linew',3);
        % indicate station name 
        text(size(STAEuse,1)*0.75,i,stnm(i,:),'fontsize',12);
    end
    evtid = strcat(num2str(evt(1,1)),num2str(evt(1,5)),num2str(evt(1,6)),...
                   num2str(evt(1,7)));
    title(strcat(evtid,{'  Normalized E component of selected regular earthquake'}));
    xlim([0 length(STAEuse(:,i))]);
    ylim([0 size(stnm,1)+1]);
    xlabel(strcat('Time*',num2str(sps)));
    box on
    hold off
    
    subplot(312)
    hold on
    for i = 1: size(stnm,1)
        % plot waveform
        plot(STANuse(:,i)/max(STANuse(:,i))+i,'color',[0.5 0.5 0.5]);
        % plot p wave arrival
        plot([atpind(i)-otind atpind(i)-otind],[-1+i 1+i],'r--','linew',1);
        text(atpind(i)-otind - 300, i+0.5, 'p/P','fontsize',10,'color','r');
        % indicate the p wave signal window
        plot([atpind(i)-otind-wlensecbf*sps, atpind(i)-otind+wlensecaf*sps-1],[i+0.5 i+0.5], ...
             'k-','linew',3);
        % plot s wave arrival
        plot([atsind(i)-otind atsind(i)-otind],[-1+i 1+i],'b--','linew',1);
        text(atsind(i)-otind - 300, i+0.5, 's/S','fontsize',10,'color','b');
        % indicate the s wave signal window
        plot([atsind(i)-otind-wlensecbf*sps, atsind(i)-otind+wlensecaf*sps-1],[i+0.5 i+0.5], ...
             'k-','linew',3);
        % indicate station name 
        text(size(STANuse,1)*0.75,i,stnm(i,:),'fontsize',12);
    end
    title('Normalized N component');
    xlim([0 length(STANuse(:,i))]);
    ylim([0 size(stnm,1)+1]);
    xlabel(strcat('Time*',num2str(sps)));
    box on
    hold off
    
    subplot(313)
    hold on
    for i = 1: size(stnm,1)
        % plot waveform
        plot(STAUuse(:,i)/max(STAUuse(:,i))+i,'color',[0.5 0.5 0.5]);
        % plot p wave arrival
        plot([atpind(i)-otind atpind(i)-otind],[-1+i 1+i],'r--','linew',1);
        text(atpind(i)-otind - 300, i+0.5, 'p/P','fontsize',10,'color','r');
        % indicate the p wave signal window
        plot([atpind(i)-otind-wlensecbf*sps, atpind(i)-otind+wlensecaf*sps-1],[i+0.5 i+0.5], ...
             'k-','linew',3);
        % plot s wave arrival
        plot([atsind(i)-otind atsind(i)-otind],[-1+i 1+i],'b--','linew',1);
        text(atsind(i)-otind - 300, i+0.5, 's/S','fontsize',10,'color','b');
        % indicate the s wave signal window
        plot([atsind(i)-otind-wlensecbf*sps, atsind(i)-otind+wlensecaf*sps-1],[i+0.5 i+0.5], ...
             'k-','linew',3);
        % indicate station name 
        text(size(STAUuse,1)*0.75,i,stnm(i,:),'fontsize',12);
    end
    title('Normalized U component');
    xlim([0 length(STAUuse(:,i))]);
    ylim([0 size(stnm,1)+1]);
    xlabel(strcat('Time*',num2str(sps)));
    box on
    hold off
    

end

% keyboard