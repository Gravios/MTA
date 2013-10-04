classdef MTADxyz < MTAData
    properties 
        path        % path to file containing object
        type        % TimeSeries,TimePeriods,TimePoints,SpatialBins
        ext
        sampleRate  % Sampling rate of data
    end
    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADxyz(path,data,type,ext,sampleRate)
            Data = Data@MTAData(path,data,type,ext,sampleRate);
        end
        function Data = create(Data,varargin)
            [Session,overwrite] = DefaultArgs(varargin,{{},{}});
            
        end

        function Data = updatePath(Data,path)
            Data.path = path;
        end
        function Data = filter(Data,win)
        end
        function Data = resample(data,newSampleRate)
        end
        function Data = embed(Data,win,overlap)
        end

        function Data = subsref(Data,S)
            ni = numel(S);
            if strcmp(S(1).type,'()')&&ni==1,
                if numel(S.subs)==1,
                    Data = Data.data(S.subs);
                elseif size(S.subs,1)==1,
                    Data = Data.data(S.subs{1},S.subs{2},S.subs{3});
                else
                    Data = SelectPeriods(Data.data(:,S.subs{2},S.subs{3}),S.subs{1},'c');
                end
                return
            end
            for n = 1:ni,
                if isa(Data,'MTADxyz'),
                    Data = builtin('subsref',Data,S(n));
                else
                    Data = subsref(Data,S(n));
                end
            end
        end
    end
    
end