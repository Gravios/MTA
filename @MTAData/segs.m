function varargout = segs(Data,varargin)
%function data = segs(Data,start_points,segment_length,if_not_complete)
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
[start_points,segment_length,if_not_complete] = ...
    DefaultArgs(varargin,{1:Data.size(1),round(.5*Data.sampleRate),nan});
oriDataSize = Data.size(2:5);
if numel(nargout)>1,
    varargout = cell(1,2);
else
    varargout = cell(1,1);
end
[varargout{:}] = reshape(GetSegs(Data.data,start_points,segment_length,if_not_complete),[segment_length,numel(start_points),oriDataSize]);

end
