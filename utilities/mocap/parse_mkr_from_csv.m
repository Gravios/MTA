function parse_mkr_from_csv(meta)
% 
% COMPLETE SOON gosh darn it 
%
% ASSUMES the csv labels are contiguous and follow the following convention
% <quaternion> <quaternion> <quaternion> <quaternion> <position> <position> <position>
% 
% REQUIRES meta.subjects


rb.name = 'RatW';
rb.nMarkers = 4;

% COLLECT subjects and their rigidbody components

rboNames = af(@(s) {s.rb.alias}, meta.subjects);
rboNames = cat(2,rboNames{:});

for take = 1:numel(meta.motive.csv),
% EXTRACT marker data from csv
    fPosition = fullfile(path.processed.ephys,meta.motive.csv{take});
    fid = fopen(fPosition,'r');

% GET header
    hdr = fgetl(fid);
    hdr = regexp(hdr,',','split');

% GET total frames 
% GET total exported frames 
    numFrames = str2double(hdr{find(~cellfun(@isempty,regexp(hdr,'Total Frames in Take')))+1});
    numExportedFrames = str2double(hdr{find(~cellfun(@isempty,regexp(hdr,'Total Exported Frames')))+1});

% SKIP empty line after header
    fgetl(fid);

% GET column headers
    clabels    = fgetl(fid);
    modelNames = fgetl(fid); modelNames = regexp(modelNames,',','split');
    elementId  = fgetl(fid); elementId  = regexp(elementId, ',','split');
    dataType   = fgetl(fid); dataType   = regexp(dataType,  ',','split');
    dimId      = fgetl(fid); dimId      = regexp(dimId,     ',','split');

    
    
    

    rigidBodyDataSize = 45;
    targetDataSize = 54;

    ndims = 3;
    numBodyMarkers = 3;

    markerData = nan([numFrames, numBodyMarkers.*ndims]);
    rigidBodyData = nan([numFrames, rigidBodyDataSize]);

    tic;
    for f = 1:numFrames,
        tdata = [fgetl(fid),','];
        if length(tdata)<10,
            continue;
        end
        cind = 0;
        ccind = [];
        ccount = 1;
        while ccount <= rigidBodyDataSize,
            cind = cind + 1;                                 
            ccount = ccount + double(tdata(cind)==',');
            if tdata(cind)==',',
                ccind(end+1) = cind;
            end
        end
        %% SOMETHING WRONG HERE (I may have fixed this... not sure)
        rowRigidBodyData = cellfun(@str2double,regexp(tdata(1:cind),',','split'));
        rigidBodyData(f,:) = rowRigidBodyData(1:end-1);
        rowMarkerData = cellfun(@str2double,regexp(regexprep(regexprep(tdata((cind-1+find(tdata(cind:end)~=',',1,'first')):end),'(,,)+',','),'(,,)+',','),',','split'));
        markerData(f,1:end+numel(rowMarkerData)-1-size(markerData,2)) = rowMarkerData(1:end-1);
        if mod(f,10000)==0,
            disp(f);
            toc;
        end
    end


fclose(fid);

frames = rigidBodyData(:,1);
timestamps = rigidBodyData(:,2);
rigidBodyData(:,1:2) = [];

save(fullfile(path.processed.xyz,meta.mazeName,[meta.sessionName,'.pos.csv_preproc.mat']),...
     'frames','timestamps','rigidBodyData','markerData','-v7.3');


    data = reshape(data,[],numel(meta.subjects.rb),8)

% REORDER from XZY to XYZ
    data(:,:,[3,4,6,7]) = data(:,:,[4,3,7,6]);
% RESCALE from meters to milimeters
    data(:,:,5:8) = data(:,:,5:8)*1e3;
    
    subjects.name = 
    
    save(fullfile(path.processed.xyz,[meta.sessionName,'.pos.csv_rbo.mat'])
end


                        
                        
                        


