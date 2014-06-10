function [accg,tbin] = autoccg(Session,varargin)
% function [accg,tbin] = autoccg(Session,varargin)
% [units,states,binSize,halfBins,normalization] = DefaultArgs(varargin,{[],[],16,60,'count'});
[units,states,binSize,halfBins,normalization,spkopt] = DefaultArgs(varargin,{[],[],16,60,'count',[]});

if ~isa(Session,'MTASession'),
    Session = MTASession(Session);
end

Session.spk.create(Session,Session.sampleRate,states,units,spkopt);

if isempty(units)
    units = 1:size(Session.spk.map);
end

accg = zeros(1+2*halfBins,size(Session.spk.map,1));
for i = units,
    uRes = Session.spk(i);
    [tccg,tbin] = CCG(uRes,i,binSize,halfBins,Session.sampleRate,i,normalization);
    accg(:,i) = tccg;
end

