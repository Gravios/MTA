classdef MTAData < hgsetget

    properties (Abstract)
        path        % path to file containing object
        data        % data
        type        % TimeSeries,TimePeriods,TimePoints,SpatialBins
        ext
        sampleRate  % Sampling rate of data

        %xyz 
        %    path - '/data/homes/gravio/data/analysis/jg05-20120317/jg05-20120317.cof.all.xyz.mat'
        %    data - xyz(TimeSeries,marker,xyzdim)
        %    type - TimeSeries
        %    sampleRate - ~119.881035 sample/sec

        %ang 
        %    data - ang(TimeSeries,marker,marker,sphdim)
        %    type - TimeSeries
        %    sampleRate - ~119.881035 sample/sec

        %Pfs
        %    data - rateMap{unit}(xbins,ybins)
        %    type - SpatialBins
        %    sampleRate - 1/diff(bin1(1:2)) sample/mm

        %Bhv.States
        %Fet.Features
        %lfp
        %ccg
        %ufr
        %Clusters.res
        %Ripples

        %Clusters
        %    data       - ResClu(timepoints,clusterId)
        %    type       - TimePoints
        %    sampleRate - 1250,32000 or 32552 sample/sec
        

    end

    methods
        function Data = MTAData(path,data,type,ext,sampleRate)
            Data.path = path;
            Data.data = data;
            Data.type = type;
            Data.ext = ext;
            Data.sampleRate = sampleRate;
        end
        function out = save(Data,varargin)
            out = false;
            data = Data.data;
            if isempty(varargin),
                data = Data.data;
                save(Data.path,'data','-v7.3');
                out = true;
            else
                data = Data.data;
                save(varargin{1},'data','-v7.3');
                out = true;
            end
        end
        function Data = load(Data,varargin)
            if isempty(varargin),
                load(Data.path)
                Data.data = data;
            else
                Data.path = varargin{1};
                load(Data.path)
                Data.data = data;
            end
        end
        function sdim = size(Data,varargin)
            if ~isempty(varargin),
                sdim = size(Data.data,cell2mat(varargin));
            else
                sdim = size(Data.data);
            end
        end

    end

    methods (Abstract)        
        Data = create(Data,varargin);
        Data = updatePath(Data,path);
        Data = filter(Data,win);
        Data = resample(data,newSampleRate);
        Data = embed(Data,win,overlap);

    end

end