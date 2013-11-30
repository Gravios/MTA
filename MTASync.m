classdef MTASync < hgsetget

    properties 
        filename
        path        % path to file containing object
        ext
        data        % data
        origin
    end
    
    methods
        function Sync = MTASync(varargin)
            [path,filename,data,ext] = DefaultArgs(varargin,{[],[],0,'sync'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename= [filename '.' ext '.mat'];
                end
            end
            Sync.filename = filename;
            Sync.path = path;
            Sync.data = data;
            Sync.origin = data(1);
            Sync.ext = ext;
        end
        function out = save(Sync,varargin)
            [overwrite,newpath] = DefaultArgs(varargin,{0,[]});
            out = false;
            if ~exist(Sync.fpath,'file')||overwrite,                
                save(Sync.path,'Sync','-v7.3');
                out = true;
            else
                save(newpath,'Sync','-v7.3');
                out = true;
            end
        end
        function Sync = load(Sync,varargin)
            load(Sync.fpath);
        end
        
        function per = periods(Sync,varargin)
            [sampleRate,ind] = DefaultArgs(varargin,{1,':'});
            per = round((Sync.data-Sync.origin).*sampleRate);
            per(per<=0)=1;
            per = per(ind,:);
        end   
        function Sync = resync(Sync,Data,varargin)
        % Sync - MTASync
        % epochs - new periods corresponding to the xyz object
            [epochs] = DefaultArgs(varargin,{[]});
            if ~isempty(epochs)
                Sync.data = epochs./Data.sampleRate+Sync.origin;            
            end
            switch Data.type
                case 'TimeSeries',                    
                    dind = false(Data.size(1),1);
                    periods = Sync.periods(Data.sampleRate);
                    for i = 1:size(Sync.periods,1),
                        dind(periods(i,1):periods(i,2))=true;
                    end 
                    Data.data(dind==0,:,:,:,:) = 0;
                    Data.data = Data(periods(1):periods(end),:,:,:,:);
                case 'TimePeriods'
                    Data.data = round(IntersectRanges(Data.data,...
                                                Sync.periods(Data.sampleRate))...
                                                -(Sync.data(1)-Sync.origin)*Data.sampleRate);
                    Data.data(Data.data==0) = 1;
                case 'TimePoints'
                    %Data.data = SelectPeriods
            end
        end
        
        function DataCopy = copy(Data)
        % Make a copy of a handle object.
        % Instantiate new object of the same class.
            DataCopy = feval(class(Data),[]);
            % Copy all non-hidden properties.
            p = properties(Data);
            for i = 1:length(p)
                DataCopy.(p{i}) = Data.(p{i});
            end
        end
        function fpath = fpath(Sync)
            fpath = fullfile(Sync.path,Sync.filename);
        end
        function Sync = updatePath(Sync,path)
            Sync.path = path;
        end
        function Sync = updateFilename(Sync,filename)
            Sync.filename = filename;
        end
        
        function Sync = subsref(Sync,S)
            ni = numel(S);
            if strcmp(S(1).type,'()')&&ni==1,
                if numel(S.subs)==0,
                    Sync = Sync.data;
                else
                    Sync = builtin('subsref',Sync.data,S);
                end
                return
            end
            for n = 1:ni,
                if isa(Sync,'MTASync'),
                    Sync = builtin('subsref',Sync,S(n:end));
                    break
                else
                    subsref(Sync,S(n:end));
                end
            end
        end        
    end

end