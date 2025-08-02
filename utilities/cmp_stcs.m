function labelStats = cmp_stcs(Trial,Stc1,Stc2,varargin)
[states,displayStcs] = DefaultArgs(varargin,{{'walk','rear','turn','pause','groom','sit'},false},true);

% $$$ Trial = 'jg05-20120317.cof.all';
% $$$ Trial = 'Ed03-20140624.cof.all';
% $$$ Trial = 'Ed01-20140707.cof.all';
% $$$ states = {'walk','rear','turn','pause','groom','sit'};
% $$$ %Stc1 = 'nn0317_PP';
% $$$ Stc1 = 'nn0624Ed_PP';
% $$$ Stc2 = 'nn0624jg_PP';
% $$$ 
% $$$ Stc1 = 'nn0707Ed_PP';
% $$$ Stc2 = 'nn0624Ed_PP';
% $$$ 
% $$$ Stc1 = 'hand_labeled_rev1_Ed';
% $$$ Stc2 = 'hand_labeled_rev1_jg';
% $$$ 
% $$$ Stc1 = 'NN0317';
% $$$ Stc2 = 'hand_labeled_rev1_jg';
% $$$ Stc1 = 'MTAC_BATCH-+seh+fet_mis_SR_12_NORM_1_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_100_NI_100_NN_multiPN_RAND_WSBNT-wrnpmsa'
% $$$ 
% $$$ 
% $$$ Stc1 = 'NN0317_PP';
% $$$ Stc2 = 'hand_labeled_rev1_Ed';
% $$$ 
% $$$ Stc1 = 'hand_labeled_rev2_Ed';
% $$$ Stc2 = 'hand_labeled_rev2_jg';


Trial = MTATrial.validate(Trial);

if ischar(Stc1)
    Stc1 = Trial.load('stc',Stc1);
end

if ischar(Stc2)
    Stc2 = Trial.load('stc',Stc2);
end

xyz = Trial.load('xyz');

shl = MTADxyz('data',double(0<stc2mat(Stc1,xyz,states)),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(0<stc2mat(Stc2,xyz,states)),'sampleRate',xyz.sampleRate); 

aper = Trial.stc{'a'};
aper.cast('TimeSeries');
aper.resample(xyz);

eper = Trial.stc{'e'};
eper.cast('TimeSeries');
eper.resample(xyz);

ind = any(shl.data,2)&any(ysm.data,2)&logical(aper.data)&~logical(eper.data);

tcm = confmat(shl(ind,:),ysm(ind,:)); % #DEP: netlab
labelStats.confusionMatrix = round(tcm./xyz.sampleRate,2);
labelStats.precision = [round(diag(tcm)./sum(tcm,2),4).*100]';
labelStats.sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
labelStats.accuracy = sum(diag(tcm))/sum(tcm(:));

if displayStcs
    hfig = figure(30230230);
    sp(1) = subplot(211);
    imagesc([1:size(shl,1)]./shl.sampleRate,1:numel(states),shl.data')
    sp(2) = subplot(212);
    imagesc([1:size(ysm,1)]./ysm.sampleRate,1:numel(states),ysm.data')       
    
    set(sp,'TickDir','out');
    set(sp,'YTick',0:numel(states));
    set(sp,'YTickLabels',states);

    linkaxes(sp,'xy');
end
