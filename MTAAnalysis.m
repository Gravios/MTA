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

    methods 
        function Anal = MTAAnalysis(parameters)
        end
    end
end