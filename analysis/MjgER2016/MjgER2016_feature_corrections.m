% req20211104
%     Tags: features correction 
%     Status: Active
%     Type: Analysis
%     Author: Justin Graboski
%     Final_Forms: NA
%     Project: General
%     Description: more empirical method for correcting deviant features


% phzCorrection
tid = 29;
Trial = Trials{tid};
unitSet = units{tid};
unitSet = unitsInts{tid};
lfp = Trial.load('lfp',sessionList(tid).thetaRefGeneral);
phz = lfp.phase([5,12]);
spk = Trial.load('spk',lfp.sampleRate,'theta-groom-sit',unitSet);
spk = Trial.load('spk',lfp.sampleRate,'rear&theta-groom-sit',unitSet);

figure();
rose(phz(spk(121)),36);

for u = 1:numel(unitSet)
    mphz(u) = circ_mean(phz(spk(unitSet(u))));
end
figure();
rose(mphz+pi/4,16);



%hbangCorrection
for tid = 1:30
Trial = Trials{tid};
hbang = fet_hbp_hba(Trial);
roll = fet_roll(Trial);
xyz = preproc_xyz(Trial,'trb');
ang = create(MTADang,Trial,filter(copy(xyz),'ButFilter',4,1.5,'low'));
vxy = vel(filter(copy(xyz),'ButFilter',4,1.5,'low'),'hcom',[1,2]);
hvfl = fet_href_HXY(Trial,[],[],'trb',[],0);
vang = circ_dist(circshift(ang(:,'hcom','nose',1),-1),circshift(ang(:,'hcom','nose',1),1)).*xyz.sampleRate;
% $$$ figure,
% $$$ histogram(vang(xper.data,1),linspace(-10,10,100),'EdgeColor','none');

xper = resample(Trial.stc{'x'}.cast('TimeSeries'),xyz);
xper.data = xper.data&sqrt(sum(sq(xyz(:,'hcom',[1,2]).^2),2))<300&abs(vang)<2;

theta = circ_mean(atan2(hvfl(xper,2)./sqrt(sum(hvfl(xper,:).^2,2)),...
           hvfl(xper,1)./sqrt(sum(hvfl(xper,:).^2,2))));
empirHrotCorrection(tid) = theta;

figDir = fullfile(MTA_FIGURES_PATH,'empircal_feature_corrections');
create_directory(fullfile(figDir));
figName = [Trial.filebase,'-hrot_hbang_roll'];

hfig = set_figure_layout(figure(1),'A4','landscape');
hfig.Units = 'normalized';
clf(hfig);
subplot2(2,3,1,1);
    hist2(hvfl(xper,:),...
          linspace(-50,100,50),...
          linspace(-80,80,50));
    caxis([0,prctile(nonzeros(get(findobj(gcf(),'Type','Image'),'CData')),99)]);
    Lines([],0,'r');
    Lines(0,[],'r');
    title('Head Speed Fwd Lat');
subplot2(2,3,2,1);
    hist2(multiprod(hvfl(xper,:),[cos(theta),-sin(theta);sin(theta),cos(theta)],[2],[1,2]),...
         linspace(-50,100,50),...
         linspace(-80,80,50))
    Lines([],0,'r');
    Lines(0,[],'r');
title(num2str(empirHrotCorrection(tid),'Hrot C %f'));
subplot2(2,3,1,2);
    set(histogram(hbang(xper,2),linspace(-pi,pi,100)),'EdgeColor','none');
    Lines(median(hbang(xper,2)),[],'r');
    empirHbangCorrection(tid) = -median(hbang(xper,2));
    title('Head-Body angle');    
subplot2(2,3,2,2);
    set(histogram(hbang(xper,2)+empirHbangCorrection(tid),linspace(-pi,pi,100)),'EdgeColor','none');
    Lines(0,[],'r');
    title(num2str(empirHbangCorrection(tid),'Hbang C %f'));
subplot2(2,3,1,3);
    set(histogram(roll(xper,1),linspace(-pi/2,pi/2,100)),'EdgeColor','none');
    Lines(median(roll(xper,1)),[],'r');
    empirRollCorrection(tid) = -median(roll(xper,1));
    title('Head Roll');
subplot2(2,3,2,3);
    set(histogram(roll(xper,1)+empirRollCorrection(tid),linspace(-pi/2,pi/2,100)),'EdgeColor','none');
    Lines(0,[],'r');
    title(num2str(empirRollCorrection(tid),'Roll C %f'));
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));    
end


% jg05
subject = 'jg05';
subjectInd = ~cellfun(@isempty,regexp({sessionList.sessionName},subject))
hfig = set_figure_layout(figure(1),'A4','landscape');
hfig.Units = 'normalized';
clf(hfig);
subplot(131);
hist(empirHrotCorrection(subjectInd),20)
    mrc = mean(empirHrotCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,'Mean Hrot C %f'));
subplot(132);
hist(empirHbangCorrection(subjectInd),20)
    mrc = mean(empirHbangCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,'Mean Hbang C %f'));
subplot(133);
    hist(empirRollCorrection(subjectInd),20)
    mrc = mean(empirRollCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,'Mean Roll C %f'));
figName = [subject,'-hrot_hbang_roll_summary'];
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));    

% if no significant changes in correction are noted then use mean for each subject
% jg05 corrections
% head direction    0.263843536305317
% head-body angle  -0.233780357672081
% roll             -0.364894511244453


% er01 corrections
subject = 'er01';
subjectInd = ~cellfun(@isempty,regexp({sessionList.sessionName},subject))
hfig = set_figure_layout(figure(1),'A4','landscape');
hfig.Units = 'normalized';
clf(hfig);
subplot(131);
hist(empirHrotCorrection(subjectInd),20)
    mrc = mean(empirHrotCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Hrot C %f']));
subplot(132);
hist(empirHbangCorrection(subjectInd),20)
    mrc = mean(empirHbangCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Hbang C %f']));
subplot(133);
    hist(empirRollCorrection(subjectInd),20)
    mrc = mean(empirRollCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Roll C %f']));
figName = [subject,'-hrot_hbang_roll_summary'];
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));    

% roll after er01-20110119 is different (marker replacement)


% er01-20110119 corrections
% head direction    0.224539353048083
% head-body angle  -0.0839221543387804
% roll             -0.0117475201093919

% er01 other
% head direction    0.251138358170965
% head-body angle  -0.0839221543387804
% roll             -0.375907022964544




% ER06 corrections
subject = 'ER06';
subjectInd = ~cellfun(@isempty,regexp({sessionList.sessionName},subject))
hfig = set_figure_layout(figure(1),'A4','landscape');
hfig.Units = 'normalized';
clf(hfig);
subplot(131);
hist(empirHrotCorrection(subjectInd),20)
    mrc = mean(empirHrotCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Hrot C %f']));
subplot(132);
hist(empirHbangCorrection(subjectInd),20)
    mrc = mean(empirHbangCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Hbang C %f']));
subplot(133);
    hist(empirRollCorrection(subjectInd),20)
    mrc = mean(empirRollCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Roll C %f']));
figName = [subject,'-hrot_hbang_roll_summary'];
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));    

% ER06 corrections
% head direction   -0.155125753662453
% head-body angle   0.15547759421534
% roll             -0.163464852468445



% Ed10 corrections
subject = 'Ed10';
subjectInd = ~cellfun(@isempty,regexp({sessionList.sessionName},subject))
hfig = set_figure_layout(figure(1),'A4','landscape');
hfig.Units = 'normalized';
clf(hfig);
subplot(131);
hist(empirHrotCorrection(subjectInd),20)
    mrc = mean(empirHrotCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Hrot C %f']));
subplot(132);
hist(empirHbangCorrection(subjectInd),20)
    mrc = mean(empirHbangCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Hbang C %f']));
subplot(133);
    hist(empirRollCorrection(subjectInd),20)
    mrc = mean(empirRollCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Roll C %f']));
figName = [subject,'-hrot_hbang_roll_summary'];
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));    

% Ed10 corrections
% head direction    0.0337460932039883
% head-body angle  -0.0769751470023147
% roll             -0.186558565753508



% jg04 corrections
subject = 'jg04';
subjectInd = ~cellfun(@isempty,regexp({sessionList.sessionName},subject))
hfig = set_figure_layout(figure(1),'A4','landscape');
hfig.Units = 'normalized';
clf(hfig);
subplot(131);
hist(empirHrotCorrection(subjectInd),20)
    mrc = mean(empirHrotCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Hrot C %f']));
subplot(132);
hist(empirHbangCorrection(subjectInd),20)
    mrc = mean(empirHbangCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Hbang C %f']));
subplot(133);
    hist(empirRollCorrection(subjectInd),20)
    mrc = mean(empirRollCorrection(subjectInd));
    Lines(mrc,[],'r');
    title(num2str(mrc,[subject,' Mean Roll C %f']));
figName = [subject,'-hrot_hbang_roll_summary'];
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));    

% jg04 corrections
% head direction   -0.136311593814964
% head-body angle   0.15984592719997
% roll              0.194087173796434


% FS03 corrections
% head direction    0
% head-body angle   0.115
% roll             -0.127


