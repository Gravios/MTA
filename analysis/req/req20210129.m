
% struct Session  
%--------------------------------------------------------------
%|  sessionName    |  String        |   'jg05-20120309'       |
%|  subjects       |  CellStr       |   {{'jg05'}}            |
%|  mazeName       |  String        |   'cof'                 |
%|  trialName      |  String        |   'all'                 |
%|  dLoggers       |  CellStr       |   {{'nlx','vicon'}}     |
%|  dPaths         |  struct        |   SEE dPaths bellow     |
%|  xyzSampleRate  |  Double        |   119.881035            |
%|  hostServer     |  String        |   'lmu'                 |
%|  dataServer     |  String        |   'lmu',                |
%|  project        |  String        |   'general'             |
%|  TTLValue       |  CONTEXT       |   SEE TTLValue bellow   |
%|  includeSyncInd |  Array[double] |   []                    |
%|  offsets        |  Array[double] |   [15,-15]              |
%|  xOffSet        |  Double        |   0                     |
%|  yOffSet        |  Double        |   0                     |
%|  stcMode        |  String        |   'default'             |
%|  rotation       |  Double        |   0                     |
%|  thetaRef       |  Array[double] |   [8:8:64]              |
%| thetaRefGeneral |  Double        |   8                     |
%--------------------------------------------------------------
%
% dPaths structure:
%----------------------------------------------------------------------------------------------
%|     FIELD       |      TYPE      |          examlple                                       |
%----------------------------------------------------------------------------------------------
%|  xyz            |  String        |   '/storage/<user>/data/processed/xyz/<subjectId>/'     |
%|  ephys          |  String        |   '/storage/<user>/data/processed/nlx/<subjectId>/'     |
%|  video          |  String        |   '/storage/<user>/data/processed/video/<subjectId>/'   |
%----------------------------------------------------------------------------------------------


meta.sessionName = 'FS03-20201222';
meta.subjects  = {'FS03'};
meta.mazeName = 'cof';
meta.trialName = 'all';
meta.dLoggers = {'WHT','CSV'};
meta.dPaths.xyz   = fullfile('/storage/gravio/data/processed/xyz/',  meta.subjects{1}{1});
meta.dPaths.ephys = fullfile('/storage/gravio/data/processed/ephys/',meta.subjects{1}{1});
meta.xyzSampleRate = 180.00;% ???
meta.hostServer = 'lmu';
meta.dataServer = 'lmu';
meta.project    = 'general';
meta.TTLValue   = [];
meta.includeSyncInd = [];
meta.offsets  = [0,0];
meta.xOffSet  = 0;
meta.yOffSet  = 0;
meta.stcMode  = 'default';
meta.rotation = 0;
meta.thetaRef = [1:11:64];
meta.thetaRefGeneral = 1;

link_session( meta.sessionName, meta.dPaths)

MTASession(meta.sessionName,                         ...
           meta.mazeName,                            ...
           true,                                        ...
           meta.TTLValue,                            ...
           meta.dLoggers,                            ...
           meta.xyzSampleRate);

path.raw.ephys = fullfile('/storage/gravio/data/raw/ephys/',meta.subjects{1},meta.sessionName);
path.processed.xyz = fullfile('/storage/gravio/data/processed/xyz/',meta.subjects{1},meta.sessionName);
fpath = dir(path.raw);fpath(1:2) = [];

%system(['tail -c +9 ',fullfile(fpath.folder,fpath.name),' > ',fullfile(fpath.folder,[filebase,'.dat'])]);
%fpath = dir(path.raw);fpath(1:2) = [];
%!cp /storage/gravio/data/project/general/jg05-20120312/jg05-20120312.xml /storage/gravio/data/raw/ephys/FS07/FS07-20201214/FS07-20201214.xml


%% PARSE rigidbody marker data ---------------------------------------------------------------------
fPosition = fullfile(path.processed.xyz,[filebase,'.pos.csv']);

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

% SET rigidbody marker indices 
targetMarkers = 31:45;
numMarkers = 5;

dstruct = substruct('()',{[1,2,31:45]});
data = nan([numFrames, numel(targetMarkers)+2]);
if numFrames == numExportedFrames,
    for f = 1:numFrames,
        data(f,:) = subsref(cellfun(@str2double,regexp(fgetl(fid),',','split')),dstruct);
    end
end

fclose(fid);


frames = data(:,1);
timestamps = data(:,2);
data(:,1:2) = [];

figure,plot3(data(:,1),data(:,3),data(:,2),'.');

% RESHAPE data to MTA standard [ time, marker,dim ]
data = permute(reshape(data,[],3,5),[1,3,2]);
data = data(:,:,[1,3,2]);

figure,plot3(data(:,1,1),data(:,1,2),data(:,1,3),'.');

% CONVERT meters -> milimeters
data = data*1000;

figure,plot3(data(:,1,1),data(:,1,2),data(:,1,3),'.');
figure,plot3(data(100,:,1),data(100,:,2),data(100,:,3),'-+');


save(fullfile(path.processed.xyz,meta.mazeName,[meta.sessionName,'.pos.data.mat']),...
     'frames','timestamps','data','-v7.3');

%% END PARSE rigidbody marker data -----------------------------------------------------------------








nind = ~isnan(data);

rowElementId  = elementId(nind);  % [hash id] ???
rowDataType   = dataType(nind);   % [rotation, position, err]
rowModelNames = modelNames(nind); % model and marker names
rowDimId      = dimId(nind);      % Dimention label [X,Y,Z,W,X,Y] oh geez guys.
rowData       = data(nind);       % the data

clear('rigidBodyObjects');
rigidBodyObjects = {struct('elementId','"none"',...
                           'modelName','none',...
                           'rowNumber',0,     ...
                           'timestamp',0,     ...
                           'frame',    0,     ...
                           'data',     []     ...
                           )                  ...
                   };
crb = 1;


% rbo structure
%  Model.name = 'RatJ';
%  Model.data = [ Rot, Rot, Rot, Rot, Pos, Pos, Pos, perMarkerError ];
%  Model.markers {
%    Marker.name = 'RatJ_1'; % RatJ_N where N = number of markers
%    Marker.data = [ Pos, Pos, Pos, quality ];
%  




data(end+1,:) = cellfun(@str2double,regexp(fgetl(fid),',','split'));




for f = 1:1%???
frame = data(1);
ftime = data(2);    
for d = 1:sum(nind),

% Check current rbo.element with the rowElementId
if strcmp(rigidBodyObjects{crb}.elementId,rowElementId(d))
    r
rbo.timestamp = ftime    
rigidBodyObjects{1}.rowNumber = f;
rigidBodyObjects{1}.elementId = rowElementId{d};
rigidBodyObjects{1}.elementId = rowElementId{t};

else
% Find existing rbo or create new
rboInd = cellfun(@(rbo,row) strcmp(rbo.elementId, row),...
                 rigidBodyObjects,...
                 repmat(rowElementId(d),size(rigidBodyObjects))...
                );
end





% $$$ while true,
% $$$     try,
% $$$         data(end+1,:) = cellfun(@str2double,regexp(fgetl(fid),',','split'));
% $$$     catch err,
% $$$         disp(err);
% $$$         break;
% $$$     end
% $$$ end



31:45



s = MTASession.validate('FS03-20201222.cof.all');
xyz = s.load('xyz');


figure,
plot(xyz.data(:,1,1)+55,xyz.data(:,1,2)-125,'.');
xlim([-500,500]);
ylim([-500,500]);
Lines([],0,'k');
Lines(0,[],'k');
Lines([],-400,'k');
Lines(-400,[],'k');
Lines([],400,'k');
Lines(400,[],'k');

figure,
plot(xyz.data(:,1,3));
xlim([-500,500]);
ylim([-500,500]);
Lines([],0,'k');
Lines(0,[],'k');
Lines([],-400,'k');
Lines(-400,[],'k');
Lines([],400,'k');
Lines(400,[],'k');

xyz.data = xyz.data+repmat(permute([55,-125,-650],[1,3,2]),[size(xyz,1),size(xyz,2),1]);


figure,
plot(xyz.data(:,1,1),xyz.data(:,1,2),'.');
xlim([-500,500]);
ylim([-500,500]);
Lines([],0,'k');
Lines(0,[],'k');
Lines([],-400,'k');
Lines(-400,[],'k');
Lines([],400,'k');
Lines(400,[],'k');


xyz.save();


label_theta(s,s.stc,1,1);

Trial = QuickTrialSetup(meta);

Trial = MTATrial.validate([meta.sessionName,'.',meta.mazeName,'.',meta.trialName]);

pft = pfs_2d_theta(Trial,'pfsArgsOverride',struct('numIter',1,'halfsample',false));

pft = pfs_2d_theta(Trial,'pfsArgsOverride',struct('numIter',1,'halfsample',false));


figure();
for p = Trial.spk.map(:,1)',
    plot(pft,p,1,'colorbar');
    title(num2str(p));
    waitforbuttonpress();
end


% NEXT step
% determine if csv file contains body marker trajectories


%% PARSE rigidbody marker data ---------------------------------------------------------------------
fPosition = fullfile(path.processed.xyz,[meta.sessionName,'.pos.csv']);

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

% $$$ data =  [];
% $$$ tdata(end+1,:) = cellfun(@str2double,regexp(fgetl(fid),',','split'));
% $$$ sdata = fgetl(fid);
% $$$ tic,cellfun(@str2double,regexp(sdata,',','split'));toc
% $$$ tic,cellfun(@str2double,regexp(sdata(1:(numel(sdata)-find(fliplr(sdata)~=',',1,'first'))),',','split'));;toc
% $$$ tic,numel(cellfun(@str2double,regexp(regexprep(fgetl(fid),',,',''),',','split'))),toc
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
    cind = 1;
    ccind = [];
    ccount = 1;
    while ccount <= rigidBodyDataSize,
        ccount = ccount + double(tdata(cind)==',');
        cind = cind + 1;                         
        if tdata(cind)==',',
            ccind(end+1) = cind;
        end
    end
    %% SOMETHING WRONG HERE
    rigidBodyData(f,:) = cellfun(@str2double,regexp(tdata(1:cind-3),',','split'));
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


%% END PARSE rigidbody marker data -----------------------------------------------------------------


ds = load(fullfile(path.processed.xyz,meta.mazeName,[meta.sessionName,'.pos.csv_preproc.mat']));

% LOOP through frames
%     HANDLE following:
%             Expected change in distance   -> use a marker-trajectoryId of previous frame as the 
%                                              marker-trajectoryId of the current frame marker, 
%                                              which minimizes the distance between frames
%             Unexpected change in distance -> generate new marker-trajectoryId
%
%     ASSIGN expected distance change in milimeters-per-second
% LOOP end
% LOOP

ds.markerData = ds.markerData * 1000;
ds.markerData = permute(reshape(ds.markerData,size(ds.markerData,1),3,9),[1,3,2]);
ds.markerData(ds.markerData==0) = nan;


figure();
plot(sqrt(sum((ds.markerData(:,1,:)-circshift(ds.markerData(:,1,:),-1)).^2,3)))

ind = 91582:91590;
sq(ds.markerData(ind,:,1))

figure();plot(sqrt(sum((ds.markerData(ind,1,:)-circshift(ds.markerData(ind,1,:),-1)).^2,3)))

trajId = nan([size(ds.markerData,1),size(ds.markerData,2)]);
xyzSampleRate= 180;
trajThresh = 10;%mm/frame
curTrajId = 1;
mdist = nan([size(ds.markerData,1),size(ds.markerData,2)]);
trajSegInds = GetSegs(-5:size(ds.markerData,1),11,1:size(ds.markerData,1),nan);
for m = 1:size(ds.markerData,2),
    f = 1;
    for f = 
        mdist(f,m) = sum(diff(sq(ds.markerData([f-1,f],m,:))).^2);        
        if mdist > trajThresh        
            curTrajId = curTrajId + 1;            
        end
        if ~isnan(mdist),
            trajId(f,m) = curTrajId;
        end
    end
end



for m = 1:9,
mdist = sum(sq(ds.markerData(f,m,:)-ds.markerData(f,1,:)).^2)
end


