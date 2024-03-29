function parse_rbo_from_csv(meta)
% function parse_rbo_from_csv(meta)
% 
% ARGIN:
%     meta - struct: contain session information including subject information
%                       necessary for this parser
% 
% - NOTICE -----------------------------------------------------------------------------------------
%
% ASSUMES the csv labels are contiguous and follow the following convention
% <quaternion> <quaternion> <quaternion> <quaternion> <position> <position> <position>
% 
% ASSUMES that the same subjects are present in all csv files for the session
%
% REQUIRES !empty meta.subjects
% REQUIRES !empty meta.csv
% 
%
% TODO ---------------------------------------------------------
%     GENERATE test cases multiple subjects with multiple RBOs
% TODO ---------------------------------------------------------
%     IF meta.csv is empty
%     GET csv list from data/raw/[meta.sessionBase]/[meta.sessionName]
%
%

% COLLECT subjects and their rigidbody components

rboAlias = af(@(s) {s.rb.alias}, meta.subjects);
rboAlias = cat(2,rboAlias{:});

cind = find(~cellfun(@isempty,meta.csv));

HEADER_OFFSET = 7;

for c = 1:numel(meta.csv),
    disp(['[INFO] MTA:utilities:mocap:parse_rbo_from_csv: Processing ', ...
          'CSV ' num2str(c), ' of ', num2str(numel(meta.csv))]);
% EXTRACT marker data from csv
    csv = meta.csv{c};
    if isempty(csv),
        continue;
    end
    
    fPosition = fullfile(meta.path.raw.xyz,csv);
    fid = fopen(fPosition,'r');

% GET header
    hdr = fgetl(fid);
    hdr = regexp(hdr,',','split');
    


% GET total frames 
% GET total exported frames 
% GET sampleRate
    numFrames = str2double(hdr{find(~cellfun(@isempty,regexp(hdr,'Total Frames in Take')))+1});
    numExportedFrames = str2double(hdr{find(~cellfun(@isempty,regexp(hdr,'Total Exported Frames')))+1});
    sampleRate = str2double(hdr{find(~cellfun(@isempty,regexp(hdr,'Export Frame Rate')))+1});    
    
% SKIP empty line after header
    fgetl(fid);

% GET column headers
    clabels    = fgetl(fid); clabels    = regexp(clabels,   ',','split');
    modelNames = fgetl(fid); modelNames = regexp(modelNames,',','split');
    elementId  = fgetl(fid); elementId  = regexp(elementId, ',','split');
    dataType   = fgetl(fid); dataType   = regexp(dataType,  ',','split');
    dimId      = fgetl(fid); dimId      = regexp(dimId,     ',','split');

    fclose(fid);
% SET rigidbody marker indices for all subjects' RBOs
% $$$     dstruct = substruct('()',{[1,2,find(~cellfun(@isempty,regexp(modelNames,['(^',strjoin(rboAlias,'$)|(^'),'$)'])))]});
% $$$ 
% $$$ 
% $$$     data = nan([numFrames, 2+8*numel(rboAlias)]);
% $$$ 
% $$$     if numFrames == numExportedFrames,
% $$$         for f = 1:numFrames,
% $$$             % NOTE - For each frame, the frame, time, and all subjects are put in the second dim of data
% $$$             %        the data matrix is later reshaped into a 3 dimension tensor <time,rbo,space>         
% $$$             data(f,:) = subsref(cellfun(@str2double,regexp(fgetl(fid),',','split')),dstruct);
% $$$             if mod(f,1000)==0
% $$$                 disp([num2str(f),' of ' num2str(numExportedFrames)])
% $$$             end
% $$$         end
% $$$     end
% $$$ 
% $$$     fclose(fid);
    numColumns = numel(find(~cellfun(@isempty,regexp(clabels,'^Rigid Body'))))+1;
    numRows = numExportedFrames+HEADER_OFFSET-1;
    csvdata = csvread(fullfile(meta.path.raw.xyz,csv),HEADER_OFFSET,0,[HEADER_OFFSET,0,numRows,numColumns]);
    
    if numFrames ~= numExportedFrames
        data = zeros([size(numFrames),57]);
        for dim = 1:size(data,2)
            data(:,dim) = interp1(csvdata(:,1),csvdata(:,dim),1:numFrames);
        end
    else
        data = csvdata;
    end    
    
% GET rbo order
    % I know, a little too complex.
    rboIndex = [];
    for s = 1:numel(rboAlias),
        rind = find(~cellfun(@isempty,regexp(modelNames,['^',rboAlias{s},'$'])),1,'first');
        assert(~isempty(rind),'MTA:utilities:mocap:parse_rbo_from_csv:RboNotFoundInCSV',rboAlias{s});
        rboIndex(s) = rind;
    end
    [~,rboIndex] = sort(rboIndex);

% SPLIT data 
    frames = data(:,1);
    timestamps = data(:,2);
    
    data(:,1:2) = [];
    
    dataRbo = data(:,find(~cellfun(@isempty,regexp(clabels,'^Rigid Body$')))-2);
    dataRbo = permute(reshape(dataRbo,[size(dataRbo,1),8,numel(rboAlias)]),[1,3,2]);
    

% REORDER subjects
    dataRbo = dataRbo(:,rboIndex,:);
% REORDER from XZY to XYZ
    dataRbo(:,:,[3,4,6,7]) = dataRbo(:,:,[4,3,7,6]);
% RESCALE from meters to milimeters
    dataRbo(:,:,5:8) = dataRbo(:,:,5:8)*1e3;
% REORDER from <Q>XYZ to XYZ<Q>
    dataRbo(:,:,:) = dataRbo(:,:,[5:8,1:4]);

% ASSIGN rbo(s) to subject(s) with correct order 
    for s = 1:numel(meta.subjects)
        subjects(s).name = meta.subjects(s).name;
        subjects(s).type = 'rbo';
        subjects(s).sampleRate = sampleRate;
        subjects(s).csv = csv;
        subjects(s).rboLabels = {};
        for r = 1:numel(meta.subjects(s).rb)
            subjects(s).rboLabels{r} = meta.subjects(s).rb(r).name;
            subjects(s).data(:,r,:) = dataRbo(:,cellfun(@strcmp,repmat({meta.subjects(s).rb(r).alias},size(rboAlias)),rboAlias),:);
        end
    end

% SAVE 
    save(fullfile(meta.path.processed.xyz,[meta.sessionName,num2str(find(c==cind),'.take_%04.f.mat')]),...
         'csv', 'frames', 'sampleRate', 'timestamps', 'numFrames', 'numExportedFrames', 'subjects', '-v7.3');

    clear('subjects'); 
    
end



                        
                        
                        


