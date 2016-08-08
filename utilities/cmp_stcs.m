function cmp_stcs(Trial,Stc1,Stc2,varargin)
[states] = DefaultArgs(varargin,{{'walk','rear','turn','pause','groom','sit'}},true)

Trial = 'jg05-20120317.cof.all';
Trial = 'Ed03-20140624.cof.all';
Trial = 'Ed01-20140707.cof.all';
states = {'walk','rear','turn','pause','groom','sit'};
%Stc1 = 'nn0317_PP';
Stc1 = 'nn0624Ed_PP';
Stc2 = 'nn0624jg_PP';

Stc1 = 'nn0707Ed_PP';
Stc2 = 'nn0624Ed_PP';

Stc1 = 'hand_labeled_rev1_Ed';
Stc2 = 'hand_labeled_rev1_jg';

Stc1 = 'hand_labeled_rev2_Ed';
Stc2 = 'hand_labeled_rev2_jg';


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
ind = any(shl.data,2)&any(ysm.data,2)&logical(aper.data);

tcm = confmat(shl(ind,:),ysm(ind,:)); % #DEP: netlab
confusionMatrix = round(tcm./xyz.sampleRate,2);
precision = round(diag(tcm)./sum(tcm,2),4).*100;
sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
accuracy = sum(diag(tcm))/sum(tcm(:));
