function varargout = segs(Data,varargin)
%function data = segs(Data,startPoints,segmentLength,ifNotComplete)
% Wraper for - > [Segs Complete] = GetSegs(x, StartPoints, SegLen, IfNotComplete);
%
% extracts multiple segments from a time series x.  x may be a vector
% or a matrix, in which case time is on the first dimension, and channel
% is the second.
%
% StartPoints give the starting times of the segments.
% SegLen gives the length of the segments.  All segments must be
% the same length, and a rectangular array is returned
% (dim 1:time within segment, dim 2: segment number, dim 3:channel);
%
% IfNotComplete specifies what to do if the start or endpoints are outside
% of the range of x.  If a value is specified, any out-of-range points
% are given that number.  If [] is specified, incomplete segments are not
% returned.  A list of complete segments extracted is returned in Complete.
% Default value for IfNotComplete: NaN.
%
[startPoints,segmentLength,ifNotComplete] = ...
    DefaultArgs(varargin,{[],round(.5*Data.sampleRate),nan});

if isempty(startPoints),
    varargout = cell([1,2]);
    return
end

oriDataSize = Data.size(2:5);
if numel(nargout)>1,
    varargout = cell(1,2);
else
    varargout = cell(1,1);
end
[varargout{:}] = reshape(GetSegs(Data.data,startPoints,segmentLength,ifNotComplete),[segmentLength,numel(startPoints),oriDataSize]);

end
