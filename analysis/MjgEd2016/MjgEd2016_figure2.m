% MjgEd2016 Figure 2 Select the behavioral feature space



Trial = MTATrial('jg05-20120317');



% SUBPLOT Section A Skeleton examples  ---------------------------------------------------------------------

slist = {'rear','sit','groom','pause','walk','turn'};

for s = 1:numel(slist)
    open(fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',['sk_',slist{s},'.fig']))

    hfig = figure(gcf);
    set(hfig,'PaperPositionMode','auto');
    set(hfig,'units','centimeters');
    set(hfig,'Position',[0,0,8,8])
    set(gca,'units','centimeters');
    set(gca,'Position',[4-4/2,4-4/2,4.5,3.5])
    set(findall(gcf,'-property','FontSize'),'FontSize',8)

    print(hfig,'-depsc2','-loose',...
          fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
                   ['fig3_A_',num2str(s),'_sk_',slist{s},'.eps']))
end


% SUBPLOT Section B t-SNE of state subset ---------------------------------------------------------------------

Trial = MTATrial.validate(Trial);
mfilename = 'req20160310_8_genOptfigs';

acc_ori = [];
acc_opt = [];

sen_ori = [];
sen_opt = [];

pre_ori = [];
pre_opt = [];

sbind = [];

dsd = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));
bs = load(fullfile(Trial.spath,'req20160310_5_genfigs.mat'));
for s = 1:5
       
    oind = [repmat([1:59],1,2)',zeros([118,1])];
    aind = oind(:,1);
    for sh = 1:117,
        oind = [oind;[circshift(aind,-sh),aind]];
    end
    slind = oind(dsd.fetInds{s},:);
    ofet =reshape(slind,[],1);
    best_inds = histc(ofet,1:59);
    [~,sbind(:,s)] = sort(best_inds,'descend');
    
    ori = load(fullfile(Trial.spath,['req20160310_4_accumStats',num2str(s),'.mat']));    
    opt = load(fullfile(Trial.spath,['req20160310_7_accumOptStats',num2str(s),'.mat']));

    acc_ori(:,s) =  cell2mat({ori.accum_acc});
    acc_opt(:,s) =  cell2mat({opt.accum_acc});    

    sen_ori(:,s) =  cell2mat({ori.accum_sen});
    sen_opt(:,s) =  cell2mat({opt.accum_sen});    

    pre_ori(:,s) =  cell2mat({ori.accum_pre});
    pre_opt(:,s) =  cell2mat({opt.accum_pre});    

end


figure
set(hfig,'PaperPositionMode','auto');
set(gcf,'units','centimeters');
set(gcf,'Position',[0,0,25,5.5])
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0,0,25,5.5])

for s = 1:5;
subplot2(1,5,1,s);
hold on,
%acc
hs = plot(acc_ori(:,s).*100);
hs.Color = [1,0.75,0.75];
hs = plot(acc_opt(:,s).*100);
hs.Color = [1,0,0];

% $$$ %sen
% $$$ hs = plot(sen_ori(:,s));
% $$$ hs.Color = [0.75,0.75,1];
% $$$ hs = plot(sen_opt(:,s));
% $$$ hs.Color = [0,0,1];
% $$$ 
% $$$ %pre
% $$$ hs = plot(pre_ori(:,s));
% $$$ hs.Color = [0.75,1,0.75];
% $$$ hs = plot(pre_opt(:,s));
% $$$ hs.Color = [0,1,0];
title(dsd.stateOrd{s})
yl = ylim;
ylim([yl(1),100]);
xlim([1,30])
set(gca,'units','centimeters');
cpos = get(gca,'Position');
set(gca,'Position',[cpos(1),cpos(2),2.5,3.5])
end
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
                     'fig3_B_accuracy.eps'))








% SUBPLOT Section C JPDF of best features--------------------------------------------------------------

dsd = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));
bs = load(fullfile(Trial.spath,'req20160310_5_genfigs.mat'));
load(fullfile(Trial.spath,'req20160310_8_genOptfigs.mat'));

% Rear Selection
% Preprocessing
s = 1;
cstate = Trial.stc{'r+s+m+p+w+n&a'};
sts = {'rear','s+m+p+w+n&a'};
stc = 'rc';
hedgs = {linspace(-2,4.5,60),linspace(-2.5,2.5,60)};
edgs  = hedgs;
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
% Preformating
hfig = figure(201603108);clf;
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,8,8])
% Plot
hist2([dsd.afet(cstate,bfets{s}(1)),dsd.afet(cstate,bfets{s}(3))],hedgs{:});
caxis([0,40])
hold on,
for i = 1:numel(sts),
    b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(1)),dsd.afet(Trial.stc{sts{i}},bfets{s}(3))];
    o = hist2(b,hedgs{:});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[5,5],'linewidth',2.5,'Color',stc(i))
end
% Labels
xlabel('Normalized Spine Pitch (A.U.)')
ylabel('Normalized Head Speed Z-axis (A.U.)')
title('Rear Vs All Excluding: Rear')
% Formating
set(gca,'units','centimeters');
set(gca,'Position',[4-3.5/2,4-3.5/2,3.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',8)
% Print
print(hfig,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
                     'fig3_C_jpdf_rear.eps'))



% Sit
% Preprocessing
s = 2;
cstate = Trial.stc{'s+m+p+w+n&a'};
sts = {'s','m+p+w+n&a'};
stc = 'rc';
hedgs = {linspace(-2.5,3,60),linspace(-2.2,2,60)};
edgs  = hedgs;
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
% Preformating
hfig = figure(201603108);clf;
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,8,8])
% Plot
hist2([dsd.afet(cstate,bfets{s}(1)),dsd.afet(cstate,bfets{s}(2))],hedgs{:});
caxis([0,40])
hold on,
for i = 1:numel(sts),
    b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(1)),dsd.afet(Trial.stc{sts{i}},bfets{s}(2))];
    o = hist2(b,hedgs{:});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[5,5],'linewidth',2.5,'Color',stc(i))
end
% Labels
xlabel('Normalized Height of Lower Body (A.U.)')
ylabel('Normalized Height of Pelvis (A.U.)')
title('Sit Vs All Excluding: Rear,Sit')
% Formating
set(gca,'units','centimeters');
set(gca,'Position',[4-3.5/2,4-3.5/2,3.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',8)
% Print
print(hfig,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
                     'fig3_C_jpdf_sit.eps'))



% Groom
% Preprocessing
s = 3;
cstate = Trial.stc{'m+p+w+n&a'};
sts = {'m','p+w+n&a'};
stc = 'rc';
hedgs = {linspace(-3,4,60),linspace(-1.5,4,60)};
edgs  = hedgs;
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
% Preformating
hfig = figure(201603108);clf;
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,8,8])
% Plot
hist2([dsd.afet(cstate,bfets{s}(1)),dsd.afet(cstate,bfets{s}(2))],hedgs{:});
caxis([0,40])
hold on,
for i = 1:numel(sts),
    b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(1)),dsd.afet(Trial.stc{sts{i}},bfets{s}(2))];
    o = hist2(b,hedgs{:});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
end
% Labels
xlabel('Normalized Height of Lower Body (A.U.)')
ylabel('Normalized Spine Curvature in XY Plane (A.U.)')
title('Groom Vs All excluding: Rear,Sit,Groom')
% Formating
set(gca,'units','centimeters');
set(gca,'Position',[4-3.5/2,4-3.5/2,3.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',8)
% Print
print(hfig,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
                     'fig3_C_jpdf_groom.eps'))

% Pause
% Preprocessing
s = 4;
cstate = Trial.stc{'p+w+n&a'};
sts = {'p','w+n'};
stc = 'rc';
hedgs = {linspace(-2,2.5,60),linspace(-2,3.5,60)};
edgs  = hedgs;
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
% Preformating
hfig = figure(201603108);clf;
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,8,8])
% Plot
hist2([dsd.afet(cstate,bfets{s}(2)),dsd.afet(cstate,bfets{s}(5))],hedgs{:});
caxis([0,40])
hold on,
for i = 1:numel(sts),
    b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(2)),dsd.afet(Trial.stc{sts{i}},bfets{s}(5))];
    o = hist2(b,hedgs{:});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
end
% Labels
xlabel('Normalized Lower Body Height (A.U.)')
ylabel('Normalized Lower Body Speed (A.U.)')
title('Pause Vs All excluding: Rear,Sit,Groom,Pause')
% Formating
set(gca,'units','centimeters');
set(gca,'Position',[4-3.5/2,4-3.5/2,3.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',8)
% Print
print(hfig,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
                     'fig3_C_jpdf_pause.eps'))


s = 5;
sts = {'walk','turn'};
stc = 'rc';
hedgs = {linspace(-1,2.5,60),linspace(-1.5,3.5,60)};
edgs  = hedgs;
[edgs{:}] = get_histBinCenters(edgs);
[X,Y] = meshgrid(edgs{:});
% Preformating
hfig = figure(201603108);clf;
set(hfig,'PaperPositionMode','auto');
set(hfig,'units','centimeters');
set(hfig,'Position',[0,0,8,8])
% Plot
hist2([dsd.afet(Trial.stc{'w+n'},bfets{s}(2)),dsd.afet(Trial.stc{'w+n'},bfets{s}(3))],hedgs{:});
caxis([0,40])
hold on,
for i = 1:numel(sts),
    b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(2)),dsd.afet(Trial.stc{sts{i}},bfets{s}(3))];
    o = hist2(b,hedgs{:});
    F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
    o = conv2(o,F,'same');
    contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
end
% Labels
xlabel('Normalized Head Speed (A.U.)')
ylabel('Normalized All Trajectories PPC (A.U.)')
title('Walk Vs Turn')
% Formating
set(gca,'units','centimeters');
set(gca,'Position',[4-3.5/2,4-3.5/2,3.5,3.5])
set(findall(gcf,'-property','FontSize'),'FontSize',8)
% Print
print(hfig,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
                     'fig3_C_jpdf_walk_turn.eps'))



% SUBPLOT Section D

slist = {'rear','sit','groom','pause','walk&turn'};

for s = 1:numel(slist)
    open(fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',['mta_tsne-',num2str(s),'.fig']))

    hfig = figure(gcf);
    set(hfig,'PaperPositionMode','auto');
    set(hfig,'units','centimeters');
    set(hfig,'Position',[0,0,8,8])
    set(gca,'units','centimeters');
    legend('location','SouthOutside')
    set(gca,'Position',[4-3.5/2,4-3.5/2,3.5,3.5])
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'box','on')    
    set(findall(gcf,'-property','FontSize'),'FontSize',8)

    print(hfig,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
                                  ['fig3_D_tnse_',slist{s},'.eps']))
end








% BEHAVIOR subtyping - Grooming -----------------------------------------------------------------------

groomPeriodsOri = [];

sessionList = get_session_list('hand_labeled');
numSessions = numel(sessionList);
referenceSessionIndex = 1;
RefTrial = MTATrial.validate(sessionList(referenceSessionIndex));

Trials = arrayfun(@(Trial) MTATrial.validate(Trial), sessionList,     'UniformOutput',false);
Stc    = cf(@(Trial) Trial.load('stc'),Trials);
Fet    = cf(@(Trial,Ref) fet_raw(Trial),Trials);
cf(@(fet,ref) fet.map_to_reference

for s = 1:numSessions,

    fet = fet_raw(Trial);
    wfs = fet.segs([],embeddingWindow);
    wfs = circshift(wfs,embeddingWindow/2,2);
    wfs = MTADxyz('data',reshape(permute(wfs,[2,1,3]),size(wfs,2),[]),'sampleRate',fet.sampleRate);
    wfs.data(isnan(wfs.data(:)))=0;    
    %fet.unity([],[],[],[5,95],Stc{'m'});
    featuresGroomSubSet = cat(1,mfet,fet(Stc{'m'},:));
    groomPeriodsOri = cat(1,groomPeriodsOri,Stc{'m'}.data);
    embeddingWindow = 64;

end
[Um,Sm,Vm] = svd(wfs(Stc{'m'},:),0);
wts = (1:embeddingWindow)./fet.sampleRate;



% DISPLAY eigen vectors
hfig = figure;
hfig.Units = 'centimeters';
hfig.Position(3:4) = [30,4];
hfig.PaperPositionMode = 'auto';
for i = 1:40,
    subplot(4,10,i);imagesc(wts,1:size(fet,2),reshape(Vm(:,i),[],size(fet,2))'),
    caxis([-0.08,0.08]);
    axis xy
end




% COMPUTE timeseries score for first 10 eigenvectors
fetM = MTADxyz('data',wfs.data*Vm(:,1),'sampleRate',fet.sampleRate);
for i = 1:40,fetM.data(:,i) = wfs.data*Vm(:,i);end
fetM.sync = fet.sync.copy;
fetM.origin = fet.origin;


states = {'walk','rear','turn','pause','groom','sit'};
sclr = 'brgcmk';
figure,
sp = [];
sp(end+1)=subplot2(10,4,[1:8],[2:4]);
plot(fetM(:,[1:4]))
sp(end+1)=subplot2(10,4,[9,10],[2:4]);
plotSTC(Stc,fet.sampleRate,'text',states,sclr,[],false);
linkaxes(sp,'x');

%3,9

v = 3;
figure,
edx = linspace(-80,80,100);
ind = Stc{'m'};
h = bar(edx,histc(fetM(ind,v),edx),'histc');
h.FaceColor = 'r';h.EdgeColor = 'r';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;
hold on
ind = Stc{'a-m-r'};
h = bar(edx,histc(fetM(ind,v),edx),'histc');
h.FaceColor = 'c';h.EdgeColor = 'c';h.FaceAlpha = 0.4;h.EdgeAlpha = 0.4;


figure,
plot(fetM(Stc{'m'},1:20))

m = fetM(Stc{'m'},1:20);

mappedX = tsne(m(1:10:end,:), [], 2, 5, 80);

figure
plot(mappedX(:,1),mappedX(:,2),'.');
