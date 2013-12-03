classdef MTASync < hgsetget

    properties 
        filename
        path        % path to file containing object
        ext
        data        % data
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
            per = round((Sync.data-Sync.Data(1)).*sampleRate);
            per(per<=0)=1;
            per = per(ind,:);
        end   
        function Sync = resync(Sync,Data,varargin)
        % Sync - MTASync
        % epochs - new periods corresponding to the xyz object
            [epochs] = DefaultArgs(varargin,{[]});
            if ~isempty(epochs)
                Sync.data = epochs./Data.sampleRate+Sync.origin;            
            elseif Data.syncOrigin == Sync(1),
                return
            end
            
            switch Data.type
                case 'TimeSeries',
                    mf = matfile(Data.fpath);
                    oldEpoch = [Data.syncOrigin,Data.syncOrigin+size(Data.data,1)/Data.sampleRate];
                    tshift = sync(1)-oldEpoch(1);
                    newOrigin = Data.syncOrigin+tshift;
                    ndata = abs(round(tshift/Data.sampleRate));                                        
                    if tshift<0,
                        if isempty(Data.treatmentRecord),
                            tsIndex = 1:ndata + round((newOrigin - Data.syncPeriods(1))*Data.sampleRate);
                            tsIneg = tsIndex<0;
                            dataSize = size(mf,'data');
                            Data.data = cat(1,zeros([sum(tsIneg),dataSize(2:end)]),mf.data(tsIndex,:,:,:,:),Data.data);
                        else
                            %Feature to be added later
                        end
                    else
                        Data.data = Data.data(ndata:end,:,:,:,:);
                    end
                    
                    %update syncOrigin to account for the shift
                    Data.syncOrigin = newOrigin;
                    

                    tshift = sync(end)-oldEpoch(end);
                    ndata = abs(round(tshift/Data.sampleRate));
                    if tshift>0,                       
                        if isempty(Data.treatmentRecord),
                            tsIndex = 1:ndata + round((newOrigin - Data.syncPeriods(1))*Data.sampleRate)+Data.size(1);
                            dataSize = size(mf,'data');
                            tsIneg = tsIndex>dataSize(1);
                            Data.data = cat(1,Data.data,mf.data(tsIndex,:,:,:,:),zeros([sum(tsIneg),dataSize(2:end)]));
                        else
                            %Feature to be added later
                        end
                    else
                        Data.data = Data.data(1:end-ndata,:,:,:,:);
                    end

                    if numel(oldEpoch)>2
                        for i = 2:numel(oldEpoch)-1,
                            tshift = sync(i)-oldEpoch(i);
                            if tshift<0 && i<numel(oldEpoch)/2,
                                if isempty(Data.treatmentRecord),
                                    
                                else
                                    Feature to be added later
                                end
                            else
                                if isempty(Data.treatmentRecord)
                                    
                                else
                                    Feature to be added later
                                end
                                
                            end
                            
                        end
                    end

                case 'TimePeriods'
                    Data.data = round(IntersectRanges(Data.data,...
                                                Data.syncPeriods*Data.sampleRate)...
                                                -(Sync.data(1)-Data.syncOrigin)*Data.sampleRate);
                    Data.data(Data.data<=0) = 1;
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