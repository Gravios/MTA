
%% New ccg stuff
Trial = MTATrial('jg05-20120315',{'CluRes'},'crt1');
%Trial = MTATrial('jg05-20120309',{'CluRes'},'all');
%Trial = MTATrial('er02-20110907',{'CluRes'},'all');
%Trial.Bhv = MTABhv(Trial,'m20130311');

%%Select immobile periods
non_immobile_bhvs = {'rear','walk','turnR','turnL','groom'};
iper = [1,size(Trial.xyz,1)];
for i = 1:length(non_immobile_bhvs),
iper = SubstractRanges(iper,Trial.Bhv.getState(non_immobile_bhvs(i)).state);
end
Trial.Bhv = Trial.Bhv.addState('i','immobile',iper);
iperdur = diff(Trial.Bhv.getState('immobile').state,1,2);
% $$$ 
% $$$ 
% $$$ 
% $$$ BHVccg0315 = {};
% $$$ nStates = length(Trial.Bhv.States);
% $$$ for i = 1:nStates,
% $$$     bhvper = Trial.Bhv.States{i}.state;
% $$$     if isempty(bhvper),continue,end
% $$$     Bccg = MTAccg(Trial,...
% $$$                   Trial.Bhv.States{i}.label,...
% $$$                   ['ccg at ' Trial.Bhv.States{i}.label ' onset and offset'],...
% $$$                   {bhvper(:,1),bhvper(:,2)},...
% $$$                   {[Trial.Bhv.States{i}.label '_onset'],[Trial.Bhv.States{i}.label '_offset']});
% $$$     BHVccg0315{i} = Bccg.filter(gausswin(3));
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ %tbccg = BHVccg;
% $$$ BHVccg = BHVccg0315;
% $$$ %BHVccg=tbccg;
% $$$ 
% $$$ hfig = figure;
% $$$ unit = 1;
% $$$ k = 1:2:nStates*2;
% $$$ while 1,
% $$$ for i = 1:nStates,
% $$$     subplot(nStates,2,k(i)),
% $$$     bar(Bccg.tbin,BHVccg{i}(:,unit,1)),
% $$$     axis tight,
% $$$     subplot(nStates,2,k(i)+1),
% $$$     bar(Bccg.tbin,BHVccg{i}(:,unit,2)),
% $$$     axis tight,
% $$$     title([Trial.Bhv.States{i}.label ' ' num2str(unit)])
% $$$ end
% $$$ 
% $$$     unit = figure_controls(hfig,unit);
% $$$     if unit==-1;return,end
% $$$ end





%%Filtered behaviors 


rper = Trial.Bhv.getState('rear').state;
interrdur = [rper(2:end,1)-rper(1:end-1,2)];
grperind = find(interrdur>3*Trial.xyzSampleRate);
%rearing offset: filtered out any with rearing offsets preceding by
%3 seconds
no_fol_rear_rearoff = rper(grperind,2); 
%rearing offset: filtered out any with rearing offsets following by
%3 seconds
no_pre_rear_rearon =  rper([1;grperind+1],1); 




wper = Trial.Bhv.getState('walk').state;
interwdur = [wper(2:end,1)-wper(1:end-1,2)];
gwperind = find(interwdur>3*Trial.xyzSampleRate);
%walking offset: filtered out any with walking offsets preceding by
%3 seconds
no_fol_walk_walkoff = wper(gwperind,2); 
%walking offset: filtered out any with walking offsets following by
%3 seconds
no_pre_walk_walkon =  wper([1;gwperind+1],1); 


[~,nfroffi,npwoni] = NearestNeighbour(no_fol_rear_rearoff,no_pre_walk_walkon);


interwrdur = abs(no_pre_walk_walkon(npwoni)-no_fol_rear_rearoff(nfroffi));
bwrperind = find(interwrdur<3*Trial.xyzSampleRate);
badr = unique(nfroffi(bwrperind));
no_fol_rearwalk_rearoff = no_fol_rear_rearoff;
no_fol_rearwalk_rearoff(badr) = [];
ys_fol_rearwalk_rearoff =no_fol_rear_rearoff(badr);


[~,nfroni,npwoffi] = NearestNeighbour(no_pre_rear_rearon,no_fol_walk_walkoff);

interwrdur = abs(no_fol_walk_walkoff(npwoffi)-no_pre_rear_rearon(nfroni));
bwrperind = find(interwrdur<3*Trial.xyzSampleRate);
badr = unique(nfroni(bwrperind));
no_pre_rearwalk_rearon = no_pre_rear_rearon;
no_pre_rearwalk_rearon(badr) = [];
ys_pre_rearwalk_rearon =no_pre_rear_rearon(badr);




group_resdescrip{1} = 'no_pre_rear_rearon';
group_restrains{1} =   no_pre_rear_rearon;

group_resdescrip{2} = 'no_fol_rear_rearoff';
group_restrains{2} =   no_fol_rear_rearoff;

group_resdescrip{3} = 'no_pre_rearwalk_rearon';
group_restrains{3} =   no_pre_rearwalk_rearon;

group_resdescrip{4} = 'no_fol_rearwalk_rearoff';
group_restrains{4} =   no_fol_rearwalk_rearoff;

group_resdescrip{5} = 'ys_pre_rearwalk_rearon';
group_restrains{5} =   ys_pre_rearwalk_rearon;

group_resdescrip{6} = 'ys_fol_rearwalk_rearoff';
group_restrains{6} =   ys_fol_rearwalk_rearoff;

group_resdescrip{7} = 'no_pre_walk_walkon';
group_restrains{7} =   no_pre_walk_walkon;

group_resdescrip{8} = 'no_fol_walk_walkoff';
group_restrains{8} =   no_fol_walk_walkoff;




pfw = MTAPlaceField(Trial,[],{'walk'},1);

pfr = MTAPlaceField(Trial,[],{'rear'},1);

for unit = 1:size(pfw.cluMap,1),
    try
        pfMaxPos(unit,:) = pfw.maxRatePos{unit}(pfw.maxRateMax{unit},:);
    end
end


Bccg = MTAccg(Trial,...
              [Trial.Bhv.States{i}.label '_filtered_frwponoff'],...
              ['ccg at ' Trial.Bhv.States{i}.label ' onset and offset filtered for precceding and following rearing or walking onset/offsets'],...
               group_restrains,...
               group_resdescrip,...
               0,... do not add a random number to filename
              'abs_dist',... partition bins by absolute distance from specified centers
               4,... number of partitions
               pfMaxPos,1);
frccg = Bccg.filter(gausswin(3));


Bccg = MTAccg(Trial,...
              [Trial.Bhv.States{i}.label '_filtered_frwponoff'],...
              ['ccg at ' Trial.Bhv.States{i}.label ' onset and offset filtered for precceding and following rearing or walking onset/offsets'],...
              {ys_pre_rearwalk_rearon,ys_fol_rearwalk_rearoff},...
              {[Trial.Bhv.States{i}.label '_onset_filt'],[Trial.Bhv.States{i}.label '_offset_filt']},...
               0)
frccg = Bccg.filter(gausswin(3));



hfig = figure;
unit = 1;
nStates = length(group_resdescrip);
k = 1;
i = 1;
while 1,
set(0
    subplot(nStates,2,k(i)),
    bar(Bccg.tbin,frccg(:,unit,1)),
    axis tight,
    subplot(nStates,2,k(i)+1),
    bar(Bccg.tbin,frccg(:,unit,2)),
    axis tight,
    title([Bccg.name ' ' num2str(unit)])
    unit = figure_controls(hfig,unit);
    if unit==-1;return,end
end

%par>1
hfig = figure(213);
pfig= figure(214);
rfig= figure(215);
unit = 1;
nStates = length(group_resdescrip);
k = 1:2:nStates*2;
g = 1:2:nStates;
set(0,'CurrentFigure',hfig);
while 1,
    clf
    for i = 1:nStates/2
        subplot(nStates/2,2,k(i)),
        imagesc(Bccg.tbin,1:size(frccg,4),sq(frccg(:,unit,g(i),:))'),axis xy
        axis tight,colorbar
        title([group_resdescrip{g(i)} ' ' num2str(unit)])
        subplot(nStates/2,2,k(i)+1),
        imagesc(Bccg.tbin,1:size(frccg,4),sq(frccg(:,unit,g(i)+1,:))'),axis xy
        axis tight,colorbar
        title([group_resdescrip{g(i)+1} ' ' num2str(unit)])
    end

    set(0,'CurrentFigure',pfig);
    clf
    pfw.plot(unit)
    hold on
    for i = 1:length(Bccg.partition_boundaries{1}(unit,:))-2
        c(i) = circle(pfMaxPos(unit,1),pfMaxPos(unit,2),Bccg.partition_boundaries{1}(unit,i+1));
    end
    set(c,'color','w');
    scatter(Trial.xyz(Bccg.posind{1},7,1),Trial.xyz(Bccg.posind{1},7,2),15,'m','fill');
    xlim([-500,500]);
    ylim([-500,500]);
    clear('c');

    set(0,'CurrentFigure',rfig);
    clf
    pfr.plot(unit)
    hold on
    for i = 1:length(Bccg.partition_boundaries{1}(unit,:))-2
        c(i) = circle(pfMaxPos(unit,1),pfMaxPos(unit,2),Bccg.partition_boundaries{1}(unit,i+1));
    end
    set(c,'color','w');
    scatter(Trial.xyz(Bccg.posind{1},7,1),Trial.xyz(Bccg.posind{1},7,2),15,'m','fill');
    xlim([-500,500]);
    ylim([-500,500]);
    clear('c');

    set(0,'CurrentFigure',hfig);
    unit = figure_controls(hfig,unit);
    if unit==-1;return,end
end


%par>1
hfig = figure;
display_mode = 'display';
unit = 1;
nStates = Bccg.partitions;
k = 1:ncol:nStates*ncol;
i = 1;
ncol = 3;
while 1,
    for i = 1:nStates
        subplot(nStates,ncol,k(i)),
        bar(Bccg.tbin,frccg(:,unit,1,i)),axis xy
        axis tight,
        subplot(nStates,ncol,k(i)+1),
        bar(Bccg.tbin,frccg(:,unit,2,i)),axis xy
        axis tight,
    end
    subplot(nStates,ncol,k+2),hold on
    pfw.plot(unit)
    for i = 1:length(Bccg.partition_boundaries{1}(unit,:))-2
        c(i) = circle(pfMaxPos(unit,1),pfMaxPos(unit,2),Bccg.partition_boundaries{1}(unit,i+1));
    end
    set(c,'color','w');
    scatter(Trial.xyz(Bccg.posind{1},7,1),Trial.xyz(Bccg.posind{1},7,2),15,'m','fill');
    xlim([-500,500]);
    ylim([-500,500]);
    clear('c');
    title([Bccg.name ' ' num2str(unit)])

    %% START - Controls
    switch display_mode
      case 'report'
        %% ReportFig
        set(gcf,'CurrentCharacter','n');
        reportfig(hfig,[Trial.filebase '.ccg.rfilt4p'],[],['unit: ' num2str(unit)],[],0);
        print(hfig,'-depsc2','-r200',['/data/homes/gravio/mrep/rfilt4p/' ...
                            Trial.name '.rfilt4p.' ...
                            num2str(unit) '.eps'])
        print(hfig,'-dpng',['/data/homes/gravio/mrep/rfilt4p/' ...
                            Trial.name '.rfilt4p.' ...
                            num2str(unit) '.png'])
        unit = unit+1;
        if unit > size(pfw.cluMap,1), return,end
      case 'display'
        %% manual controls
        unit = figure_controls(unit);
        if unit==-1;return,end
    end
    %% END - Controls
end


%%quick rast
rperlfp = round((Trial.Bhv.getState('rear').state-1)./Trial.xyzSampleRate.*Trial.lfpSampleRate);
rpos = 782953;
rstart = rpos-round(3*Trial.xyzSampleRate);
rstop =  rpos+round(3*Trial.xyzSampleRate);
resind = find(Trial.res>rstart/Trial.xyzSampleRate*Trial.lfpSampleRate&Trial.res<rstop/Trial.xyzSampleRate*Trial.lfpSampleRate);
res = Trial.res(resind);
clu = Trial.clu(resind);
uclu = unique(clu);
figure,plot(res,clu,'.')
Lines((rpos-1)/Trial.xyzSampleRate*Trial.lfpSampleRate,[],'k');
figure,plot(res,clu,'.')
Lines((rpos-1)/Trial.xyzSampleRate*Trial.lfpSampleRate,[],'k');


rperlfp = round((Trial.Bhv.getState('rear').state-1)./Trial.xyzSampleRate.*Trial.lfpSampleRate);
twin =250000;
resind = find(Trial.res<twin);
res = Trial.res(resind);
clu = Trial.clu(resind);


rperlfpn = rperlfp(rperlfp(:,2)<twin,:);


figure,plot(res,clu,'.')
Lines(rperlfpn(:,1),[],'k');
Lines(rperlfpn(:,2),[],'r');


clulist = [39,51,55,63,68]

clures  = [33000]

clulist = [36,55,63,68,74];
clures  = [101200]


clulist = [16,88,93,114]
clures  = [101200]


clulist = [11,16,88,90,93]
clures  = [170200]


Lines(repmat(res,2,1),repmat(clu,2,1)+cat(1,ones(size(clu,1),1).*-.5,ones(size(clu,1),1).*.5),'k');
Lines(rperlfpn([14,15],1),[],'k');
Lines(rperlfpn([14,15],2),[],'r');

rpos = 174000/Trial.lfpSampleRate*Trial.xyzSampleRate;
rstart = rpos-round(5*Trial.xyzSampleRate);
rstop =  rpos+round(5*Trial.xyzSampleRate);
resind = find(Trial.res>rstart/Trial.xyzSampleRate*Trial.lfpSampleRate&Trial.res<rstop/Trial.xyzSampleRate*Trial.lfpSampleRate);

res = Trial.res(resind);
clu = Trial.clu(resind);

myres = [];
myclu = [];
clulist = [13,34,88,93,90,16,11]
for i = 1:length(clulist),
tmyres = res(clu==clulist(i));
tmyclu = i.*ones(length(tmyres),1);
myres = cat(1,myres,tmyres);
myclu = cat(1,myclu,tmyclu);
end




figure('color',[1,1,1]),
y = repmat(myclu,1,2)+cat(2,ones(size(myclu,1),1).*-.5,ones(size(myclu,1),1).*.5);
for i = 1:size(myres,1),
Lines(myres(i)./Trial.lfpSampleRate,y(i,:),'r');
end
Lines(rperlfpn([14,15],1)./Trial.lfpSampleRate,[],'g',[],2);
Lines(rperlfpn([14,15],2)./Trial.lfpSampleRate,[],'k',[],2);
Lines([],[4.5],'b',[],1);
ylim([0,8])
xlabel('Time (Seconds)')
ylabel('Neurons')

figure,plot(res,clu,'.')
Lines((rpos-1)/Trial.xyzSampleRate*Trial.lfpSampleRate,[],'k');





%% UFR vs ACC
if isempty(Trial.ufr),
    Trial = Trial.load_ufr;
end


win = 121;
comb = Filter0(gausswin(win)./sum(gausswin(win)),Trial.com(Trial.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'})));
comh = Filter0(gausswin(win)./sum(gausswin(win)),Trial.com(Trial.Model.rb({'head_back','head_left','head_front','head_right'})));

combv = sqrt(sum(diff(comb(:,[1,2])).^2,2));
comhv = sqrt(sum(diff(comh(:,[1,2])).^2,2));

ufrsegs = {};
chvsegs = {};
for i = 1:length(group_restrains),
ufrsegs{i} = GetSegs(Trial.ufr(:,16),round((group_restrains{i}-1)./Trial.xyzSampleRate.*Trial.lfpSampleRate)-2*Trial.lfpSampleRate,4*Trial.lfpSampleRate,0);
chvsegs{i} = GetSegs(abs(diff(comhv)),round(group_restrains{i}-2*Trial.xyzSampleRate),round(4*Trial.xyzSampleRate),0);
end
ufrslen = size(ufrsegs{1},1);
chvslen = size(chvsegs{1},1);

for i = 1:length(group_restrains),
ufrsegs{i} = ufrsegs{i}(round(Trial.lfpSampleRate/Trial.xyzSampleRate/2:Trial.lfpSampleRate/Trial.xyzSampleRate:ufrslen),:);
end


figure,plot(ufrsegs{6},chvsegs{6},'.')

for i = 1:length(group_restrains),
figure,plot(mean(ufrsegs{i},2),mean(chvsegs{i},2),'.')
title(group_resdescrip{i}),
end


for i = 1:length(group_restrains),
figure,plot(sum(ufrsegs{i},2),sum(chvsegs{i},2),'.')
title(group_resdescrip{i}),
end
