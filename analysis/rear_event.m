function rear_event(Session,varargin)
[trialName] = DefaultArgs(varargin,{'all'});

Session = MTASession(Session);
Session = MTATrial(Session,{},trialName);
rsts = {'pre_rear_onset','pre_peri_rear_onset','dia_rear_onset','post_peri_rear_onset','post_rear_onset'; ...
      'pre_rear_offset','pre_peri_rear_offset','dia_rear_offset','post_peri_rear_offset','post_rear_offset'};

for evt = 1:size(rsts,1),
for sts = 1:size(rsts,2),
MTAPlaceField(Session,[],rsts{evt,sts},1);
end
end



% $$$ 
% $$$     switch event_type
% $$$       case 'onset'
% $$$         rsts = {'pre_rear_onset','pre_peri_rear_onset','dia_rear_onset','post_peri_rear_onset','post_rear_onset'}
% $$$       case 'offset'
% $$$         rsts = {'pre_rear_offset','pre_peri_rear_offset','dia_rear_offset','post_peri_rear_offset','post_rear_offset'};
% $$$     end
