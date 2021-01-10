function [accg,tbin] = autoccg(Session,varargin)
% function [accg,tbin] = autoccg(Session,varargin)
% [units,states,binSize,halfBins,normalization] = DefaultArgs(varargin,{[],[],16,60,'count'});
[units,states,binSize,halfBins,normalization,spk] = DefaultArgs(varargin,{[],[],16,60,'count',[]});

if ~isa(Session,'MTASession'),
    Session = MTASession(Session);
end

if isempty(spk),
    spk = Session.load('spk',Session.sampleRate,states,units);
end

if isempty(units)
    units = 1:size(spk.map);
end

accg = zeros(1+2*halfBins,size(spk.map,1));
for i = units,
    uRes = spk(i);
    [tccg,tbin] = CCG(uRes,i,binSize,halfBins,Session.sampleRate,i,normalization);
    accg(:,i) = tccg;
end

