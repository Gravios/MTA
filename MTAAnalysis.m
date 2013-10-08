classdef MTAAnalysis < hgsetget
% Oh god here we go again
%

    properties (Abstract)
        path
        ext
        parameters
        mdata
        data
    end

    methods (Abstract)
        Anal = MTAAnalysis(parameters, funcHandle, varargin)
        Anal= genPath(Anal)
    end
end