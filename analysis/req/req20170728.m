
sessionList = 'hand_labeled';
Trials = af(@(Trial) MTATrial.validate(Trial), get_session_list(sessionList));
stc = cf(@(t) t.load('stc'), Trials);
%stc = cf(@(t) t.load('stc','msnnN0+hand_labeled'), Trials);

s = 1
xyz = cf(@(Trial) Trial.load('xyz'),                         Trials);
vxy = cf(@(x)     x.vel({'spine_lower','head_front'},[1,2]), xyz   );
for s = 1:numel(Trials), vxy{s}.data(vxy{s}.data<1e-3,:) = 1e-3;end
for s = 1:numel(Trials), vxy{s}.data = log10(vxy{s}.data);end
cf(@(x) x.filter('ButFilter',3,5,'low'),xyz);;
ang = cf(@(t,x) create(MTADang,t,x),Trials,xyz);

 

ind = cf(@(stc) [stc{'a-s-m'}],stc);

%ind = cf(@(s,w,t) [s.get_state_transitions(t,{'pause+walk','rear'},1,w)],...
%         stc, xyz, Trials);

% $$$ shift = [0,30];
% $$$ for s = 1:numel(Trials), 
% $$$     ind{s}(sum([ind{s}+repmat(shift,size(ind{s},1),1)<=0, ...
% $$$                 ind{s}+repmat(shift,size(ind{s},1),1)>size(xyz{s},1)],2)>0,:)=[];
% $$$ end

figure
for marker = {'spine_lower','pelvis_root','spine_middle'},
    pdb = cf(@(ang,ind) circ_dist(circshift(ang(ind,marker,'spine_upper',2),-1),...
                                  circshift(ang(ind,marker,'spine_upper',2), 1)),...
             ang,ind);

    pdh = cf(@(ang,ind) circ_dist(circshift(ang(ind,'head_back','head_front',2),-1),...
                                  circshift(ang(ind,'head_back','head_front',2), 1)),...
             ang,ind);

    adb = cf(@(ang,ind) circ_dist(circshift(ang(ind,marker,'spine_upper',1),-1),...
                                  circshift(ang(ind,marker,'spine_upper',1), 1)),...
             ang,ind);

    adh = cf(@(ang,ind) circ_dist(circshift(ang(ind,'head_back','head_front',1),-1),...
                                  circshift(ang(ind,'head_back','head_front',1), 1)),...
             ang,ind);

% $$$ figure,
% $$$ plot([pdb,pdh])
% $$$ plot([adb,adh])

    p = xcorr(cat(1,pdb{:}),cat(1,pdh{:}),120,'coeff');
    a = xcorr(cat(1,adb{:}),cat(1,adh{:}),120,'coeff');
    c = xcorr(sqrt((cat(1,adb{:})+cat(1,pdb{:})).^2),sqrt((cat(1,adh{:})+cat(1,pdh{:})).^2),120,'coeff');

    subplot(131),hold on
    plot([-120:120]./xyz{1}.sampleRate,a)
    grid('on');xlim([-1,1])
    subplot(132),hold on
    plot([-120:120]./xyz{1}.sampleRate,p)
    grid('on');xlim([-1,1])
    subplot(133),hold on
    plot([-120:120]./xyz{1}.sampleRate,c)
    %plot([-120:120]./xyz{1}.sampleRate,sqrt(p.^2+a.^2))
    grid('on');xlim([-1,1])
end


phaseFrequencyRangeLow = [1.2,6];
phaseFrequencyRangeHigh = [6,12];
% LOAD features
features = cf(@(t) fet_bref(t), Trials);
featuresPhaseLow  = cf(@(f,frq) f.phase(frq), features,repmat({phaseFrequencyRangeLow},[1,numel(Trials)]));
featuresPhaseHigh = cf(@(f,frq) f.phase(frq), features,repmat({phaseFrequencyRangeHigh},[1,numel(Trials)]));

% FILTER features
ffet = cf(@(f) f.copy(), features);
cf(@(f) f.filter('ButFilter',3,[1.2,6],'bandpass'),ffet);

% FIND hip swing peaks
fetLocMins = cf(@(f) LocalMinima(-abs(f(:,17)),5,-2.5), ffet);
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

