function [STAEtmr,STANtmr,STAUtmr,STAttmr] = ...
    rdhinettmr(sacpath,date,stainfo,evtall,tmrall,sps,lo,hi,npo,npa,wlensec,pltflag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the function to read the tremor data from Japan Hinet (after or before instrument
% response removal), bandpass filtering, resample if needed. By checking and avoiding all
% regular earthquake events (from 1 min before origin time to 15 min after p arrival), 
% choose the segment that is continuous, without inference of events not at the both
% ends of the day, length as requested
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

stnm = char(stainfo{1,1});    % station name, string
stlo = stainfo{1,2};    % station longitude
stla = stainfo{1,3};    % station latitude

for i = 1: size(stnm,1)
    
    clear tmpe tmpn tmpu tmpef tmpnf tmpuf tmpefd tmpnfd tmpufd
    datestr = num2str(date);
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
    
    [tmpe,hdr,~,~,time]=readsac(staepath,0,'l');   %%% Function 'readsac'
    [tmpn,~,~,~,~]=readsac(stanpath,0,'l');
    [tmpu,~,~,~,~]=readsac(staupath,0,'l');
    
    spstrue = round(1./hdr.DELTA);

    %%% ALWAYS filter first, then resample, to avoid time-domain waveform distortion
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
    
end

dt = 1/sps;
tmrtmp = tmrall(tmrall(:,1)==date, :);    % all tremor detections in that day
tmrdetind = [];    % indexes of all tremor detections in that day
for i = 1: size(tmrtmp)
    tmrstsec = tmrtmp(i,5)*3600;
    tmredsec = (tmrtmp(i,5)+1)*3600;
    tmrstind = find(abs(tmrstsec-time)<=dt,1,'first');
    tmredind = find(abs(tmredsec-time)<=dt,1,'first')-1;
    tmrdetind = [tmrdetind tmrstind:1:tmredind];
end
tmrdetind = unique(tmrdetind);      % just in case

dumpsec = 15*60;    % secs after p arrival
dumosec = 1*60;     % secs before origin time

evttmp = evtall(evtall(:,1)==date, :);  % all regular events in that day
otsec = [];
for i = 1: size(evttmp,1)
    otsec(i) = evttmp(i,5)*3600+evttmp(i,6)*60+evttmp(i,7);     % origin time
    otind(i) = find(abs(otsec(i)-time)<=dt,1,'first');    % origin time index
end

for i = 1: size(stnm,1)
    
    dumind = [];

    for j = 1: size(evttmp,1)
        
        tt=taupTime('ak135', evttmp(j,10),'p','evt',[evttmp(j,9) evttmp(j,8)],...
                    'sta',[stla(i) stlo(i)]);
        if isempty(tt)
            fprintf('No p at station %s, using P instead. \n',stnm(i,:));
            tt=taupTime('ak135', evttmp(j,10),'P','evt',[evttmp(j,9) evttmp(j,8)],...
                        'sta',[stla(i) stlo(i)]);
        end
        ttp = tt.time;   % travel time p
        atpsec = otsec(j)+ttp;    % arrival time p
        
        dumoind = find(abs(otsec(j)-dumosec-time)<=dt,1,'first');    % start index before origin
        dumpind = find(abs(atpsec+dumpsec-time)<=dt,1,'first');     % end index after p arrival
        dum = (dumoind: 1 :dumpind)';   % these indexes are segments containing EQ records potentially
        dumind = [dumind; dum];     % dummy index to be thrown away
        
        
    end
    
    dumind = dumind';
    dumind = [1:10*60*sps dumind size(time,1)-10*60*sps+1:size(time,1)];  % also throw away 10 min at both ends of the day
    dumind = dumind';
    
    dumind = unique(dumind);    % sort without repetitions  
    
    tmrallind = setdiff(1:1:size(STAt,1), dumind);   % this is all indexes safe to tremor estimate    
    
    tmrind = intersect(tmrallind, tmrdetind,'sorted');  % the final indexes should also have tremor detections
    
%     wlensec = 400;    % length of win in sec for consecutive tremor
    wlen = wlensec*sps;     % length in samples
    wlenmax = (wlensec+wlensec)*sps;    % safest length in samples with buffers at both ends
    % find the consecutive indexes and group them
    groupind = mat2cell( tmrind, 1, diff( [0, find(diff(tmrind) ~= 1), length(tmrind)] )) ;
    
    for j = 1: size(groupind,2)
        if size(groupind{j},2) > wlenmax   % the one should have the length > tmrwmax 
            tmp = groupind{j};
            tmp = tmp';
            prefind(:,i) = tmp(round(wlensec*0.5*sps): round(wlensec*1.5*sps)-1);  % find the 1st eligible consecutive window for tremor estimate
            break
        end
    end

    % window contains tremor + noise
    STAEtmr(:,i) = STAE(prefind(:,i), i);    
    STANtmr(:,i) = STAN(prefind(:,i), i);
    STAUtmr(:,i) = STAU(prefind(:,i), i);
    STAttmr(:,i) = STAt(prefind(:,i), i);
end


if pltflag
    
    [scrsz, res] = pixelperinch(2);

    figure
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 11;   % maximum height allowed is 11 inches
    set(gcf,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    
    subplot(311)
    hold on
    for i = 1: size(stnm,1)
        plot(STAttmr(:,i),STAEtmr(:,i)/max(STAEtmr(:,i))+i,'color',[0.5 0.5 0.5]);
        text(size(STAEtmr,1)*0.75,i,stnm(i,:),'fontsize',12);
        sttime = strcat(datestr,{' '},num2str(STAttmr(1,i)),{' s'});
        text(size(STAEtmr,1)*0.1, i+0.5, sttime,'fontsize',10,'color','b');

    end
    title(strcat(datestr,{'  Normalized E component of selected tremor segment'}));   
    ylim([0 size(stnm,1)+1]);
    xlabel('Time (s)');
    box on
    hold off
    
    subplot(312)
    hold on
    for i = 1: size(stnm,1)
        plot(STAttmr(:,i),STANtmr(:,i)/max(STANtmr(:,i))+i,'color',[0.5 0.5 0.5]);
        text(size(STANtmr,1)*0.75,i,stnm(i,:),'fontsize',12);
        sttime = strcat(datestr,{' '},num2str(STAttmr(1,i)),{' s'});
        text(size(STANtmr,1)*0.1, i+0.5, sttime,'fontsize',10,'color','b');
    end
    title('Normalized N component');
    ylim([0 size(stnm,1)+1]);
    xlabel('Time (s)');
    box on
    hold off
    
    subplot(313)
    hold on
    for i = 1: size(stnm,1)
        plot(STAttmr(:,i),STAUtmr(:,i)/max(STAUtmr(:,i))+i,'color',[0.5 0.5 0.5]);
        text(size(STAUtmr,1)*0.75,i,stnm(i,:),'fontsize',12);
        sttime = strcat(datestr,{' '},num2str(STAttmr(1,i)),{' s'});
        text(size(STAUtmr,1)*0.1, i+0.5, sttime,'fontsize',10,'color','b');
    end
    title('Normalized U component');
    ylim([0 size(stnm,1)+1]);
    xlabel('Time (s)');
    box on
    hold off
    
    
    figure
    widin = 8;  % maximum width allowed is 8.5 inches
    htin = 11;   % maximum height allowed is 11 inches
    set(gcf,'Position',[scrsz(1)+1*scrsz(3)/10 scrsz(2)+scrsz(4)/20 widin*res htin*res]);
    
    hold on
    i=6;
    plot(STAE(:,i),'color',[0.5 0.5 0.5]);
    for j = 1: size(evttmp,1)
        plot([otind(j) otind(j)],[-max(STAE(:,i)) max(STAE(:,i))],'r--');
    end
    plot([prefind(1,i) prefind(1,i)],[-max(STAE(:,i)) max(STAE(:,i))],'--','color',[0 1 1]);
    plot([prefind(end,i) prefind(end,i)],[-max(STAE(:,i)) max(STAE(:,i))],'--','color',[0 1 1]);
    text(size(STAE,1)*0.75,0,stnm(i,:),'fontsize',12);
    
    title(strcat(datestr,stnm(i,:),{'  E component'}));    
    xlabel(strcat('Time*',num2str(sps)));
    box on
    hold off

    
end

% keyboard
















