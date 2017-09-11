

% LOAD Trials 
sessionList = 'hand_labeled';
Trials = af(@(Trial) MTATrial.validate(Trial), get_session_list(sessionList));
Trials(2) = [];
numTrials = numel(Trials);

% SET figure directories
OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_4/parts';
mkdir(fullfile(OwnDir,FigDir));


% LOAD State collections 
stc = cf(@(t) t.load('stc'), Trials);
%stc = cf(@(t) t.load('stc','msnnN0+hand_labeled'), Trials);

% LOAD Positions
xyz = cf(@(Trial) Trial.load('xyz'),                         Trials);

fxyz = cf(@(x) x.copy, xyz);
cf(@(x) x.filter('ButFilter',3,2.5,'low'),fxyz);;
% LOAD log10 Speeds
vxy = cf(@(x)     x.vel({'spine_lower','head_front'},[1,2]), fxyz   );
for s = 1:numel(Trials), vxy{s}.data(vxy{s}.data<1e-3) = 1e-3;end
for s = 1:numel(Trials), vxy{s}.data = log10(vxy{s}.data);end

% FILTER xyz before computing angles 
cf(@(x) x.filter('ButFilter',3,8,'low'),xyz);;

% COMPUTE intermarker angles
ang = cf(@(t,x) create(MTADang,t,x),Trials,xyz);

% COMPUTE body referenced body decomposition 
features          = cf(@(t)     fet_bref(t)      , Trials  );
featuresPhaseLow  = cf(@(f,frq) f.phase([1.2, 6]), features);
featuresPhaseHigh = cf(@(f,frq) f.phase([6  ,12]), features);
% FILTER features
ffet = cf(@(f) f.copy(), features);
cf(@(f) f.filter('ButFilter',3,[1.2,6],'bandpass'),ffet);

affet = cf(@(f) f.data, ffet);
affet = cat(1,affet{:});

paffet = cf(@(f) f.data, featuresPhaseLow);
paffet = cat(1,paffet{:});


avxy = cf(@(f) f.data, vxy);
avxy = cat(1,avxy{:});

vbins = discretize(avxy(:,1),linspace(0.5,2,100));
cm = jet(100);


steps = LocalMinima(affet(:,17),10,-5);
steps(isnan(vbins(steps)))=[];

wper = cf(@(s) cast(s{'w'},'TimeSeries')  ,stc);
wind = cf(@(w) w.data, wper);a
wind = cat(1,wind{:});
steps(~wind(steps))=[];



figure,
scatter(circ_dist(paffet(steps,17),paffet(steps,19)),...
        circ_dist(paffet(steps,17),paffet(steps,21)),...
        15,cm(vbins(steps),:));
daspect([1,1,1])




ind = repmat({':'},[1,numel(Trials)]);; %cf(@(stc) [stc{'a-s-m'}],stc);

ind = cf(@(stc) [stc{'a-s-m'}],stc);
%ind = cf(@(s,w,t) [s.get_state_transitions(t,{'pause+walk','rear'},1,w)],...
%         stc, xyz, Trials);

% $$$ shift = [0,30];
% $$$ for s = 1:numel(Trials), 
% $$$     ind{s}(sum([ind{s}+repmat(shift,size(ind{s},1),1)<=0, ...
% $$$                 ind{s}+repmat(shift,size(ind{s},1),1)>size(xyz{s},1)],2)>0,:)=[];
% $$$ end

pdh = cf(@(ang,ind) circ_dist(circshift(ang(ind,'head_back','head_front',2),-1),...
                              circshift(ang(ind,'head_back','head_front',2), 1)),...
         ang,ind);
adh = cf(@(ang,ind) circ_dist(circshift(ang(ind,'head_back','head_front',1),-1),...
                              circshift(ang(ind,'head_back','head_front',1), 1)),...
         ang,ind);

pdb = {}
adb = {}
markers = {'spine_lower','pelvis_root','spine_middle'};
for m = 1:numel(markers),
    pdb(end+1,:) = cf(@(ang,ind) circ_dist(circshift(ang(ind,markers{m},'spine_upper',2),-1),...
                                  circshift(ang(ind,markers{m},'spine_upper',2), 1)),...
                      ang,ind);
    adb(end+1,:) = cf(@(ang,ind) circ_dist(circshift(ang(ind,markers{m},'spine_upper',1),-1),...
                                  circshift(ang(ind,markers{m},'spine_upper',1), 1)),...
                      ang,ind);
end



s = 1;
tind = [277,283];
figure
sp = [];
for m = 1:numel(markers),
    sp(end+1) = subplot2(2,8,1,1:5);hold('on')
    plot([1:size(pdb{m,s},1)]./xyz{s}.sampleRate,pdb{m,s},'LineWidth',1);
    plot([1:size(pdh{s},1)]./xyz{s}.sampleRate,pdh{s},'LineWidth',1);
    sp(end+1) = subplot2(2,8,2,1:5);hold('on')
    plot([1:size(adb{s},1)]./xyz{s}.sampleRate,adb{m,s},'LineWidth',1);
    plot([1:size(adh{s},1)]./xyz{s}.sampleRate,adh{s},'LineWidth',1);
    p = xcorr(cat(1,pdb{m,:}),cat(1,pdh{:}),240,'coeff');
    a = xcorr(cat(1,adb{m,:}),cat(1,adh{:}),240,'coeff');
    subplot2(2,8,1,[7:8]); hold('on'); grid('on');
    plot([-240:240]./xyz{1}.sampleRate,a,'LineWidth',1)
    xlim([-2,2])
    subplot2(2,8,2,[7:8]); hold('on'); grid('on');
    plot([-240:240]./xyz{1}.sampleRate,p,'LineWidth',1)
    xlim([-2,2]);
end
linkaxes(sp,'x')
xlim(sp(1),tind);
xlim(sp(2),tind);
exampleTimePeriodStr = strjoin(num2str(tind),'-'
FigName = ['ccg_head_body_angular_velocity,'_',exampleTimePeriodStr];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



figure,sp = [];
sp(end+1) = subplot2(3,1,1:2,1);
hold('on');
plot(adb{1,s});
plot(adh{1,s});
sp(end+1) = subplot2(3,1,3,1);
plotSTC(Stc{s});
linkaxes(sp,'x')

figure,sp = [];
sp(end+1) = subplot2(3,1,1:2,1);
hold('on');
plot(ang{s}(:,1,4,2));
plot(ang{s}(:,4,7,2));
sp(end+1) = subplot2(3,1,3,1);
plotSTC(Stc{s});
linkaxes(sp,'x')


% ANTON 20170822 HEAD BODY PITCH 
% JPDF of head body pitch 
hbpitch = cf(@(t)  fet_HB_pitch(t),  Trials);
cf(@(f,t,r) f.map_to_reference_session(t,r),...
   hbpitch, Trials, repmat({'jg05-20120317.cof.all'},[1,numTrials]));
for s = 1:numTrials, hbpitch{s}.data(~nniz(xyz{s}),:,:) = 0;end

ind = cf(@(s)  s{'w+n+p'},  stc);
      cf(@(i)  cast(i,'TimeSeries'),  ind);
ind = cf(@(i)  i.data==1,  ind);
ind = cat(1, ind{:});

ahbpitch = cf(@(f)  f.data,  hbpitch);
ahbpitch = cat(1, ahbpitch{:});

axyz = cf(@(f)  f.data,  xyz);
axyz = cat(1, axyz{:});


hfig = figure(); 
hist2([ahbpitch(ind,1),axyz(ind,5,3)],linspace(-0.3,0.4,100),linspace(50,140,100));
caxis([0,1000])


hfig = figure(); 
hist2(ahbpitch(ind,:),linspace(-0.3,0.4,100),linspace(-1,.5,100));
caxis([0,1000])

% FIND threshold for head pitch
% $$$ figure();
% $$$ subplot2(3,1,1:2,1); plot(hbpitch{1}(:,2));
% $$$ subplot2(3,1,  3,1); plotSTC(stc{1});
% $$$ linkaxes(findobj(gcf,'Type','axes'),'x');

% SET threshold for head pitch
pThresh = repmat({-0.2},[1,numTrials]);

pind = cf(@(p,t)  p(:,2)>t,  hbpitch, pThresh);
pind = cat(1, pind{:});
pper = ThreshCross(pind,0.5,20);
pper(ind(pper(:,1))==0|ind(pper(:,2))==0,:) = [];

pSegs = GetSegs(ahbpitch,pper(:,1)-30,60);
pSegs(:,sq(abs(diff(pSegs([1,end],:,2))))<0.3,:) = [];
                    
% OVERLAY Trajectories at head pitch threshold 

xl = [-0.3,0.5];
yl = [-1,0.4];
hfig = figure(); 
subplot2(1,3,1,1);grid('on')
hist2(ahbpitch(ind,:),linspace([xl,100]),linspace([yl,100]));
caxis([0,1000])
xlabel('body pitch (rad)');ylabel('head pitch (rad)');
title('JPDF walk+turn+pause');
subplot2(1,3,1,2);grid('on')
hold('on');
cm = jet(60);
for i = 1:200
    scatter(pSegs(:,i,1),pSegs(:,i,2),10,cm,'filled');
end 
xlim(xl); ylim(yl);
xlabel('body pitch (rad)');ylabel('head pitch (rad)');
title('Half second Trajectories change in pitch > 0.3 rad');
subplot2(1,3,1,3);grid('on')
hist2(reshape(pSegs,[],2),linspace([xl,100]),linspace([yl,100]));
caxis([0,100]);
xlabel('body pitch (rad)');ylabel('head pitch (rad)');
title('JPDF of Trajectories change in pitch > 0.3 rad');
ForAllSubplots('daspect([1,1,1])');
FigName = ['jpdf_head_body_pitch'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));


% HYPOTHESIS - upward head movements occupy a different route
%              compared to downward head movements




% COORDINATION of lateral rostro-caudal axial kinematics

% FIND hip swing peaks
fetLocMins = cf(@(f) LocalMinima(-abs(f(:,17)),5,-2.5), ffet);
% SELECT steps within 
fetLocMinsWalk = cf(@(i,s) SelectPeriods(i,s{'w'},'d',1), fetLocMins,stc);



figure,
plot(featuresPhase{s}(fetLocMinsWalk{s},17),...
     featuresPhase{s}(fetLocMinsWalk{s},23),'.')

figure,
plot(featuresPhaseLow{s}(fetLocMinsWalk{s},23),...
     featuresPhaseHigh{s}(fetLocMinsWalk{s},23),'.')


eds = linspace(-1,1.5,50);
[~,abin] = histc(ang{s}(:,'head_back','head_front',2),eds);
ac = jet(50);

figure,
scatter(featuresPhase{s}(fetLocMinsWalk{s},17),...
        featuresPhase{s}(fetLocMinsWalk{s},25),10,ac(abin(fetLocMinsWalk{s}),:),'filled')



figure,hold on, for s = 1:6, plot(featuresPhase{s}(fetLocMinsWalk{s},23),ang{s}(fetLocMinsWalk{s},'head_back','head_front',2),'.b'),end
figure,hold on, for s = 1:6, plot(featuresPhase{s}(fetLocMinsWalk{s},23),xyz{s}(fetLocMinsWalk{s},'head_front',3),'.b'),end
figure,hold on, for s = 1:6, plot(featuresPhase{s}(fetLocMinsWalk{s},25),vxy{s}(fetLocMinsWalk{s},2),'.b'),end




out = {};
for s = 1:6,
out{s} = hist2([featuresPhase{s}(fetLocMinsWalk{s},25),ang{s}(fetLocMinsWalk{s},'head_back','head_front',2)],-pi:.3:pi,eds);
end
figure,imagesc(-pi:.3:pi,eds,sum(cat(3,out{:}),3)'),axis xy


out = {};
for s = 1:6,
out{s} = hist2([circ_dist(featuresPhase{s}(fetLocMinsWalk{s},23),featuresPhase{s}(fetLocMinsWalk{s},17)),vxy{s}(fetLocMinsWalk{s},2)],-pi:.2:pi,0.5:0.05:2);
end
figure,imagesc(-pi:.3:pi,0.5:0.1:2,sum(cat(3,out{:}),3)'),axis xy




out = {};
for s = 1:6,
out{s} = hist2([featuresPhase{s}(fetLocMinsWalk{s},23),xyz{s}(fetLocMinsWalk{s},'head_front',3)],-pi:.2:pi,10:5:180);
end
figure,imagesc(-pi:.2:pi,0:5:180,bsxfun(@rdivide,sum(cat(3,out{:}),3),sum(sum(cat(3,out{:}),3),1))'),
axis xy
caxis([0,.14])



edy = repmat({linspace(0,320,250)},[1,numel(Trials)]);
edx = repmat({linspace(-1,1.8,250)},[1,numel(Trials)]);

s = 3
%ind = [stc{s}{'a-m-s'}];
ind = cf(@(s) [s{'p'}],stc);
out = cf(@(a,x,i,ex,ey) hist2([a(i,'spine_middle','spine_upper',2),...
                               x(i,'head_front',3)],ex,ey),...
         ang,xyz,ind,edx,edy);

figure,
imagesc(edx{1},edy{1},sum(cat(3,out{:}),3)');
caxis([0,250])
axis xy



% rhm-ncp phase difference as func of head pitch and head height




[dstate,mdstate] = compute_cross_correlation(...
        'state',  stateTransition,...
        'sbound', sbound,...
        'overwrite',overwrite);


%% mutinfo stuff
OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_4/parts';

sbound = -240:240;
randomIndex = [];
rng(20170728);
while numel(randomIndex)<1000
    randomIndex = randi([-1e5,1e5],[1,10000]);
    randomIndex(abs(randomIndex)<240) = [];
end
randomIndex = randomIndex(1:1e3);

overwrite = false;
stateTransitions = {{'pause','walk'},{'pause','turn'},{'pause+walk','rear'}};
stateTransition = stateTransitions{1};
for stateTransition = stateTransitions,
    stateTransition = stateTransition{1};
    [dstate,mdstate] = compute_time_lagged_mutual_information(...
        'state',  stateTransition,...
        'sbound', sbound,...
        'overwrite',overwrite);
    [dRnd,mdRnd] = compute_time_lagged_mutual_information(...
        'state',  stateTransition,...
        'sbound', randomIndex,...
        'overwrite',overwrite);

    [mixy,sixy] = max(permute(dstate,[2,3,1]),[],3); 
    sixy = sq(sixy)-ceil(numel(sbound)/2);
    mthresh = [sq(nanmean(dRnd))+sq(nanstd(dRnd)).*5];
    mixy(mixy<mthresh|mixy<mthresh') = nan;
    sixy(isnan(mixy)) = nan;


    sp = gobjects(0);
    cp = gobjects(0);
    hfig = figure;
    hfig.Units = 'centimeters';


    subplot(121); 
    [sp(end+1),cp(end+1)] = imagescnan(mixy',[0,2],'linear',1,[0.5,0.5,0.5],'colorMap',@bone); 
    axis('xy');title(['Maximum mutual information']);
    Lines(5.5,[],'k');Lines(10.5,[],'k')
    Lines([],5.5,'k');Lines([],10.5,'k')

    subplot(122); 
    [sp(end+1),cp(end+1)] = imagescnan(sixy'./mdstate.sampleRate,[-0.35,0.35],'linear',1,[0.5,0.5,0.5]); 
    axis('xy');title(['Time lag of maximum mutual information']);
    Lines(5.5,[],'k');Lines(10.5,[],'k')
    Lines([],5.5,'k');Lines([],10.5,'k')

    hfig.Position(3:4) = [16,5];
    for s = 1:numel(sp), 
        sp(s).Units = 'centimeters';
        sp(s).Position(3:4) = [3,3];
        cp(s).Units = 'centimeters';
        cp(s).Position([1,4]) = [sp(s).Position(1)+3.2,3];
    end

    FigFilebase = fullfile(OwnDir,FigDir,['mis-time_lagged-new-',strjoin(stateTransition,'2')]);
    print(hfig,'-depsc2',[FigFilebase,'.eps']);
    print(hfig,'-dpng',  [FigFilebase,'.png']);
    pause(0.1);
    delete(hfig)

end

figure,
plot(ixy(:,6,8))



