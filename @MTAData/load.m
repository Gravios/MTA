function Data = load(Data,varargin)
% Data = load(Data,varargin)
% Loads data from disk into MTAData object
% 
% inputs: 
%    Session: 
%    filename
%    syncshift
%
[Session,filename,syncshift] = DefaultArgs(varargin,{[],[],0});
if ~isempty(filename),
    if exist(fullfile(Data.path,filename),'file')
        Data.filename = filename;
        ds = load(Data.fpath);                    
    else
        files = dir(Data.path);
        re = ['\.' Data.ext '\.'];
        stsFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
        
        if ~isempty(filename)
            re = ['\.' filename '\.'];                   
        else
            re = ['\.' Data.label '\.'];                   
        end
        Data.filename = stsFileList{~cellfun(@isempty,regexp(stsFileList,re))};
        ds = load(Data.fpath);
    end
else
    ds = load(Data.fpath);
end


for field = fieldnames(ds)',
    field = char(field);
    Data.(field) = ds.(field);
end


switch class(Session)
  case 'MTASession'
    Data.sync.sync = Session.sync.copy;
    Data.resync(Session);                    
  case 'MTATrial'
    Data.sync.sync = Session.sync.copy;
    Data.resync(Session);                    
  case 'double'
    if ~isempty(Session),
        mf = matfile(Data.fpath);
        d = ones(1,5);
        if Data.isempty,
            dsize(1:numel(size(mf,'data'))) = size(mf,'data');
        else
            dsize(1:numel(Data.size)) = Data.size;
        end
        d(1:numel(dsize)) = dsize;
        for i = 1:size(Session,1),
            Data.data(Session(i,1):Session(i:2),:,:,:,:) = ds.data((Session(i,1):Session(i:2))-syncshift,1:d(2),1:d(3),1:d(4),1:d(5));
        end
        if isfield(ds,'sampleRate'),
            Data.sampleRate = ds.sampleRate;
        end
    else
        Data.data = ds.data;
    end
end
end
