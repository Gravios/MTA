function xyOccupancy(Session,varargin)
%function xyOccupancy(FileBase,varargin)
%[xyzSubsetName,type,res,TrackingMarker,xyzSamplingRate] = DefaultArgs(varargin,{'','mov',400,'head_front',119.881035});
% all recommended caxis = 
% mov recommended caxis = [0 5]
% rear recommended caxis = [60 120]
[states,nbins] = DefaultArgs(varargin,{'head',400});

[stsp stateLabel] = Session.statePeriods(states);
stsxyz = SelectPeriods(Session.xyz, round((stsp-Session.syncPeriods(1,1))/Session.lfpSampleRate*Session.xyzSampleRate+1), 'c', 1);
stsxyz = reshape(stsxyz,[],size(Session.xyz,2),size(Session.xyz,3));

[o xb yb p] = hist2(sq(stsxyz(:,Session.Model.gmi(Session.trackingMarker),1:2)),nbins,nbins);
imagesc(o)
caxis([0,10])