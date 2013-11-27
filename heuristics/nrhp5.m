function [nrhp_state,nrhp_feature] = nrhp5(Session,method,varargin)

nrhp_feature = [];

rper = Session.Bhv.getState('rear').state;
hper = Session.Bhv.getState('head').state;

switch method

  case 'exclusion'
    [trim] = DefaultArgs(varargin,{5});
    %% Non-rearing head velocity thresholded periods
    trimLength = round(Session.xyzSampleRate*trim);%samples
    temp_non_rearing_periods = [rper(:,1)-trimLength,rper(:,2)+trimLength];
    tnrp = [[1;temp_non_rearing_periods(:,2)],[temp_non_rearing_periods(:,1);size(Session.xyz,1)]];
    nrp = tnrp(diff(tnrp,1,2)>0,:);
    nrhp_state = IntersectRanges(nrp,hper);

end
