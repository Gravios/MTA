function Session = load_optitrack(Session,xyzSampleRate)
% function Session = load_optitrack(Session,xyzSampleRate)
%
% process and store optitrack data which has been labeled in
% ViconIQ in an MTASession object

% LOAD the xyz data from all c3d.mat files associated with the Session
if isempty(xyzSampleRate),
    [xyzDataChunks, markers, xyzSampleRate] = concat_optitrack_files(Session);            
else
    [xyzDataChunks, markers] = concat_optitrack_files(Session);
end


% LOAD VSK to get marker names,colors and connections for creating the model
vsk_path = fullfile(Session.spath, Session.maze.name, [Session.name '-' Session.maze.name '.vsk']);
if exist(vsk_path,'file'),
    Session.model = MTAModel(vsk_path,'-vsk');
else
    warning(['VSK file associated with this session was ' ...
             'not found. \nCreating a general marker model.']);
    Session.model = MTAModel(markers{1},'-mar');
end


% CORRECT the size of each block.
% unfortunately functions within the previous steps of the
% processing pipeline randomly clip data. I added an extra point
% just in case. Now it has to be stripped out. :(
for i = 1:size(xyzDataChunks,1)
    nChunks = sum(~cellfun(@isempty,xyzDataChunks(i,:)));
    for j = 1:nChunks,
        markerInds = Session.model.gmi(markers{i,j});
        markerInds(markerInds==0) = [];
        missingMarkerInds = find(~ismember(1:Session.model.N,markerInds));        
        if ~isempty(missingMarkerInds),
            if ~isempty(markerInds),
                xyzDataChunks{i,j}(:,markerInds,:) = xyzDataChunks{i,j};
                xyzDataChunks{i,j}(:,missingMarkerInds,:) = ...
                    zeros([size(xyzDataChunks{i,j},1),numel(missingMarkerInds),3]);
            else
                xyzDataChunks{i,j} = ...
                    zeros([size(xyzDataChunks{i,j},1),Session.model.N,3]);
            end
        elseif any( markerInds~=1:Session.model.N ),
            xyzDataChunks{i,j}(:,markerInds,:) = xyzDataChunks{i,j};
        end
        
        if j<nChunks-1
            if xyzDataChunks{i,j}(end,1,1)==xyzDataChunks{i,j+1}(1,1,1)
                xyzDataChunks{i,j}(end,:,:) = [];
            end
            if j==1, continue; end
            xyzDataChunks{i,1} = cat(1,xyzDataChunks{i,1},xyzDataChunks{i,j});
        else
            % Check if final length matches expected length
            % cellfun('length',{1,:}}
            if ndims(xyzDataChunks{i,j})<3
                xyzDataChunks{i,j} = permute(xyzDataChunks{i,j},[3,1,2]);
            end
            
            xyzDataChunks{i,1} = cat(1,xyzDataChunks{i,1},xyzDataChunks{i,j});                
        end
    end
end

xyzDataChunks = xyzDataChunks(:,1)';

% Create the xyzPeriods and concatinate all trials into one xyz array
xyz = double(cell2mat(xyzDataChunks(~cellfun(@isempty,xyzDataChunks))'));

% Create XYZ data object


syncPeriods = cellfun(@length,xyzDataChunks);
syncPeriods(syncPeriods==0)=[];
syncPeriods = [ cumsum([1,syncPeriods(1:end-1)]);cumsum(syncPeriods)]'...
               ./xyzSampleRate;
syncPeriods(1) = 0;

Session.sync = MTADepoch(Session.spath,[Session.filebase '.sync.mat'],syncPeriods([1,end]),1,0,0,[],[],[],'sync');


syncPeriods = MTADepoch([],[],syncPeriods,1,Session.sync.copy,0);



Dxyz = MTADxyz(Session.spath,Session.filebase,xyz,xyzSampleRate,...
               syncPeriods,0,Session.model);
Dxyz.save;

Session.xyz = Dxyz;

Session.ang = MTADang(Session.spath,Session.filebase,[],xyzSampleRate,...
                      Session.xyz.sync,Session.xyz.origin,Session.model);

Session.fet = MTADfet(Session.spath,...
                      [],...
                      [],...
                      [],...
                      Session.sync.copy,...
                      Session.sync.data(1),...
                      []);                  

%% MTAStateCollection object holds all behavioral sets of periods
Session.stc = MTAStateCollection(Session.spath,Session.filebase,'default',[],[],1);
Session.stc.updateSync(Session.sync);
Session.stc.updateOrigin(0);

Session.save
