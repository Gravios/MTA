
%% FIG2 Features and Segmentation of behaviors-------------------|
%  A: Trajectories of behaving rat                               |
%  B: Feature matrix demonstrating the features we used          |
%  C: Labels corresponding to the hand labeled data              |
%  D: Labels corresponding to the neural network labeled data    |
%  E: t-sne dimensionality reduction method                      |
%  F: Tabel of labeling stats between and animals and between    |
%     labelers                                                   |
%________________________________________________________________|



Trial = MTATrial('jg05-20120317');
Stc = Trial.load('stc','hand_labeled_rev2_jg'); 
figPath = '/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_2';
pPad = [0,0];

%exPer = [26664,27100];
%exPer = [26000,26480];
%exPer = [25200,25580];
%exPer = [55200,60000];
exPer = [51549, 55145];
%exPer = [129470,139470];




xyz = Trial.load('xyz').filter('ButFilter',3,50);
ang = create(MTADang,Trial,xyz);
%stateColors = 'brcgym';


hfig = figure(38239385);clf
set(hfig,'position',[1016,111,775,840];
set(hfig,'paperposition',[0,0,775/100,840/100])
%% Fig:3:A Skeleton examples


axes('Position', [0.1300,0.6569,0.7750,0.350]);hold on;

% $$$ set(gca,'CameraPositionMode', 'manual'                    ,...
% $$$ 	'YLim', [-200 200],...
% $$$         'XLim', [-300 400],...
% $$$ 	'ZLim', [0 300],...
% $$$ 	'CameraPosition',     [-1909.49 3535.01 1621.08],...
% $$$ 	'CameraTargetMode',   'manual'                    ,...
% $$$ 	'CameraTarget',       [50 0 150]                ,...
% $$$  	'CameraUpVectorMode', 'manual'                    ,...
% $$$  	'CameraUpVector',     [0 0 1]                   ,...
% $$$  	'CameraViewAngleMode','manual'                  ,...
% $$$  	'CameraViewAngle',    [6.88708]);
% $$$ daspect([1,1,1])
set(gca,'CameraPositionMode', 'manual'                    ,...
	'XLim', [-300 400],...
	'YLim', [-300 400],...
	'ZLim', [0 300],...
        'CameraPosition', [2050.7 4543.64 1748.25],...
	'CameraPositionMode','manual',...
	'CameraTarget',[50 50 150],...
	'CameraTargetMode','manual',...
	'CameraUpVector',[0 0 1],...
	'CameraUpVectorMode','manual',...
	'CameraViewAngle',[6.31812],...
	'CameraViewAngleMode','manual')
daspect([1,1,1]);
hold on        
pMode = 'line';  %'surface';

plotSkeleton(Trial,xyz,exPer(1)+1500,pMode,ang,[0,500]);% rear
plotSkeleton(Trial,xyz,exPer(1)+2000,pMode,ang,[0,300]);% rear
plotSkeleton(Trial,xyz,exPer(1)+2300,pMode,ang,[0,0]);% rear


zlim([0,300]);



% Fig:2:B - feature matrix
[fet,flabels,fdisc] = fet_tsne_rev10(Trial);
axes('Position',[ 0.1300,0.4,0.7750,0.2000])
ts = (1:fet.size(1))./fet.sampleRate;
per = round(exPer./xyz.sampleRate)+pPad;
ind = ts>per(1)&ts<per(2);
ufet = unity(fet);
imc = imagesc(ts(ind),1:fet.size(2),ufet(ind,:)');
caxis([-2,2]);
xlim(round(exPer./xyz.sampleRate)+pPad)
Lines(round([exPer(1)+1500,exPer(1)+2000,exPer(1)+2300]./xyz.sampleRate),[],'k');
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',1:numel(flabels),...
        'YTickLabel',flabels);
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
set(gca,'TickDir','out');


% Fig:2:C - Expert Labels
stateLabels = {'walk','rear','turn','pause','groom','sit'};
stateColors = 'brgymc';
Stc = Trial.load('stc','hand_labeled_rev2_jg');
axes('Position',[ 0.1300,0.355,0.7750,0.0400])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+pPad)
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'EXP'});
set(gca,'XTickLabelMode','manual',...
        'XTickMode','manual',...
        'XTick',[],...
        'XTickLabel',{});
set(gca,'TickDir','out');



% Fig:2:D - NN Model Labels
stateLabels = {'walk','rear','turn','pause','groom','sit'};
stateColors = 'brgymc';
Stc = Trial.load('stc','NN_multiPN-Ed03-20140625.cof.all-RAND_wsbhand_labeled_rev1-wrnpms');
axes('Position',[ 0.1300,0.31,0.7750,0.0400])
plotSTC(Stc,1,'patch',stateLabels,stateColors);
xlim(round(exPer./xyz.sampleRate)+pPad)
set(gca,'YTickLabelMode','manual',...
        'YTickMode','manual',...
        'YTick',.5,...
        'YTickLabel',{'NN'});
set(gca,'TickDir','out');

% Fig:2:E - t-Sne

axes('Position',[ 0.1300,0.025,0.38,0.250])


featureSet = 'fet_tsne_rev5';
normalize = true;
mapToRefTrial = true;
mfilename = 'req20151203';
RefTrial = MTATrial('jg05-20120317');
if mapToRefTrial, mapping =    ['-map2_' RefTrial.filebase];else mapping    = '';end
if normalize,     normStatus = '-norm';                     else normStatus = '';end
fileLoc = fullfile(Trial.path.data,'analysis',...
                   ['req20151203-hand_labeled-fet_tsne_rev5',mapping,normStatus,'.mat']);
ds = load(fileLoc);


osts = numel(ds.states);
hfig = figure(3923924);clf
hold on;
mc = ds.csmat(ds.ind,:);
for nc = 1:osts,
    nind = all(bsxfun(@eq,ds.c(nc,:),mc),2);
    h = scatter(ds.mappedX(nind,1),ds.mappedX(nind,2),2,mc(nind,:));
    try,h.MarkerFaceColor = h.CData(1,:);end
end
legend(states);
set(gca,'XTickLabel',{})
set(gca,'YTickLabel',{})


% Fig:3:F - Labeling Accuracy 








