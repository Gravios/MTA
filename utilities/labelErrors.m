function Stc = labelErrors(Trial,Stc,varargin)

v = vel(Trial.load('xyz').filter(gtwin(.25,Trial.xyz.sampleRate)),...
        Trial.trackingMarker,[1,2]);

