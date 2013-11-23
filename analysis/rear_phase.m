
Trial = MTATrial('jg05-20120316',{'CluRes',{'lfp',72}},'all');




bhvl='rear';
rper = Trial.Bhv.getState(bhvl).state;

interrdur = [rper(2:end,1)-rper(1:end-1,2)];
grperind = find(interrdur>3*Trial.xyzSampleRate);
%bhvl offset: filtered out any with bhvl offsets preceding by 3 seconds

rdur = diff(rper,1,2);
drperind = find(rdur>1.5*Trial.xyzSampleRate);

gind = sum(repmat(grperind,1,length(drperind))'==repmat(drperind,1,length(grperind)))>0;

bhvoff = rper(grperind(gind),2);
%bhvl offset: filtered out any with bhvl offsets following by 3 seconds
bhvon =  rper([1;grperind(gind)+1],1);




Bccg = MTAccg(Trial,...
              ['rear'],...
              ['ccg at ' bhvl ' onset and offset filtered: id3,d1.5'],...
               {bhvon,bhvoff},...
               {[bhvl '_onset'],[bhvl '_offset']});
frccg = Bccg.filter(gausswin(3));



%%Need theta phase
[ThPh,ThAm,TotP] = myThetaPhase(Trial.lfp);


%plot phases precession during rearing

theta = Trial.statePeriods('theta');


rz = sq(Trial.xyz(:,7,3));

GoodPeriod{1} = IntersectRanges(theta, dotdot(bhvon,'+',[-3 3]*1250));
GoodPeriod{2} = IntersectRanges(theta, dotdot(bhvoff,'+',[-3 3]*1250));

GoodPeriod{1} = IntersectRanges(theta, dotdot(bhvon,'+',[-2 2]*1250));
GoodPeriod{2} = IntersectRanges(theta, dotdot(bhvoff,'+',[-2 2]*1250));

gpz{1} = IntersectRanges(theta, dotdot(bhvon,'+',[-1.0 1.0]*1250));
gpz{2} = IntersectRanges(theta, dotdot(bhvoff,'+',[-1.0 1.0]*1250));





%figure
k = 1;
while k~=-1,
    clf
    for s=1:2,
        myph = ThPh(SelectPeriods(Trial.res(Trial.clu==k),GoodPeriod{s},'d',1,0));
        myres = SelectPeriods(Trial.res(Trial.clu==k),GoodPeriod{s},'d',1,2);
        myphz = ThPh(SelectPeriods(Trial.res(Trial.clu==k),gpz{s},'d',1,0));
        myz = rz(round(SelectPeriods(Trial.res(Trial.clu==k),gpz{s},'d',1,0)./Trial.lfpSampleRate.*Trial.xyzSampleRate));
        subplot2(5,2,1,s);
try
        [hcnt xbin ybin] =  hist2([repmat(myres,2,1)/1.250-2000,[myph*180/pi; 360+myph*180/pi]], 100, 36);
        smimage(xbin,ybin, hcnt,5);axis xy
        title(num2str(Trial.map(k,:)));
end
        subplot2(5,2,2,s);
        plot(myres/1.250-2000, myph*180/pi,'k.','MarkerSize',5); hold on
        plot(myres/1.250-2000, 360+myph*180/pi,'k.','MarkerSize',5); hold on
        axis tight
        set(gca,'YTick',[-180:180:540]); grid on
        subplot2(5,2,3,s);
try
        [hcnt xbin ybin] =  hist2([repmat(myz,2,1),[myphz*180/pi; 360+myphz*180/pi]], 100, 36);
        smimage(xbin,ybin, hcnt,5);axis xy
        title(num2str(Trial.map(k,:)));
end
        subplot2(5,2,4,s);
        plot(myz, myphz*180/pi,'k.','MarkerSize',5); hold on
        plot(myz, 360+myphz*180/pi,'k.','MarkerSize',5); hold on
        axis tight
        set(gca,'YTick',[-180:180:540]); grid on
        subplot2(5,2,5,s);
        bar(Bccg.tbin, frccg(:,k,s)); axis tight
        xlim([-3000 3000]);
    end
    k = figure_controls(gcf,k)
end

% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ GoodPeriod{1} = IntersectRanges(RUN, dotdot(rper(:,1),'+',[-1 .5]*1250));
% $$$ GoodPeriod{2} = IntersectRanges(RUN, dotdot(rper(:,2),'+',[-2 2]*1250));
% $$$ 
% $$$ s= 1;
% $$$ k= 68;
% $$$ 
% $$$ myphc = ThPh(SelectPeriods(Trial.res(Trial.clu==gi(k)),GoodPeriod{s},'d',1,0));
% $$$ myresc = SelectPeriods(Trial.res(Trial.clu==gi(k)),GoodPeriod{s},'d',1,2);
% $$$ 
% $$$ 
% $$$ 
% $$$ gpz{1} = IntersectRanges(RUN, dotdot(rper(:,1),'+',[-1 0]*1250));
% $$$ gpz{2} = IntersectRanges(RUN, dotdot(rper(:,2),'+',[-1.5 1.5]*1250));
% $$$ 
% $$$ myphzc = ThPh(SelectPeriods(Trial.res(Trial.clu==gi(k)),gpz{s},'d',1,0));
% $$$ myzc = rz(round(SelectPeriods(Trial.res(Trial.clu==gi(k)),gpz{s},'d',1,0)./Trial.lfpSampleRate.*Trial.xyzSampleRate));
% $$$ 
% $$$ 
% $$$ myphs = GetSegs(ThPh,rper(:,1),4*1250,0);
% $$$ myphs = GetSegs(ThPh,rper(:,1)-2*1250,4*1250,0);
% $$$ for i= 1:size(myphs,1),
% $$$ pval(i) = circ_rtest(myphs(i,:));
% $$$ end
% $$$ figure
% $$$ subplot2(3,1,1,1)
% $$$ plot(circ_mean(myphs'));
% $$$ subplot2(3,1,2,1)
% $$$ plot(circ_var(myphs'));
% $$$ subplot2(3,1,3,1)
% $$$ plot(pval);
% $$$ 
% $$$ 
% $$$ myphs = GetSegs(ThPh,rper(:,2)-2*1250,4*1250,1);
% $$$ for i= 1:size(myphs,1),
% $$$ pval(i) = circ_rtest(myphs(i,:));
% $$$ end
% $$$ 
% $$$ figure
% $$$ subplot2(3,1,1,1)
% $$$ plot(circ_mean(myphs'));
% $$$ subplot2(3,1,2,1)
% $$$ plot(circ_var(myphs'));
% $$$ subplot2(3,1,3,1)
% $$$ plot(pval);
% $$$ 
% $$$ 
% $$$ mylfps = GetSegs(Trial.lfp,rper(:,2)-2*1250,4*1250,1);

