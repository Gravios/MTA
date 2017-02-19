% NAME :    req20170201
% PURPOSE : Convert MTA objects to 2d mat for vasiliki
%           katsageorgiou
%

outputPath = '/storage/gravio/ownCloud/colaboration/vasiliki/';

sessionList = get_session_list('hand_labeled');

for s = sessionList
    Trial = MTATrial.validate(s);
    xyz = Trial.load('xyz');
    markerLabels = xyz.model.ml;
    behaviorState = stc2mat(Trial.stc,xyz);
    behaviorLabels = Trial.stc.list_state_attrib;
    xyz = xyz.data;
    frameRate = Trial.xyz.sampleRate;
    save(fullfile(outputPath,[Trial.filebase '-raw_data.mat']),'xyz','markers','behaviorState','behaviorLabels','frameRate');
end




RefTrial = MTATrial.validate('jg05-20120317.cof.all');
frameRate = 12;
for s = sessionList
    Trial = MTATrial.validate(s);    
    [features,featureLabels,featureDescriptions] = fet_mis(Trial,frameRate);
    behaviorState = stc2mat(Trial.stc,features);
    behaviorLabels = Trial.stc.list_state_attrib;
    featuresRaw = features.data;
    features.map_to_reference_session(Trial,RefTrial);
    featuresMapped = features.data;
    referenceTrial = RefTrial.filebase;
    save(fullfile(outputPath,[Trial.filebase '-features-fet_mis.mat']), ...
         'featuresRaw','featuresMapped','featureLabels','featureDescriptions','behaviorState','behaviorLabels','frameRate','referenceTrial');
end
