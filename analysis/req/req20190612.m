% req20190612
%    Tags: spk clustering klustakwik
%    Status: retired
%    Type: utility
%    Author: Justin Graboski
%    Final_Forms: 
%    Project: 
%    Description: reclustering after the removal of ripples


% 1. Backup files - { *.res.?, *.fet.?, *.clu.?, *.spk.? };
% 2. load spikes  - timestamp, clusterId, waveform
% 3. remove noise - where clusterId == (1 | 0 )
% 4. save files   - { *.res.?, *.spk.? };

% TODO : start ndm_pca automatically
% TODO : start KlustaKwik automatically
% TODO : rewrite backups to do incremental backups of only noise groups

%sessionName = 'g10-20130425';
sessionName = 'jg05-20120312';
sessionPath = fullfile('/storage/gravio/data/project/spkclean/',sessionName);
sessionPathOri = fullfile(sessionPath,['ori_',datestr(now,30)]);
create_directory(sessionPathOri);
bkpFileList = {'res','fet','clu','spk'};
Par = LoadPar(fullfile(sessionPath,[sessionName,'.xml']));

spikeGroupInds = 3:11;
spikeGroupInds = 4:11;



if isempty(spikeGroupInds),
% DESIGNATE all spike groups     
    spikeGroupInds = 1:numel(Par.SpkGrps);
end

for spkg = spikeGroupInds
    spkgs = num2str(spkg);
    nChannels = numel(Par.SpkGrps(spkg).Channels);
    nSamples = Par.SpkGrps(spkg).nSamples;

    
% SAVE backup of all files
    try,
        cf(@(ext) ...
           movefile(fullfile(sessionPath,[sessionName,'.',ext,'.',spkgs]),...
                    fullfile(sessionPathOri,[sessionName,'.',ext,'.',spkgs])),...
           bkpFileList);
    end
    
% LOAD spike timestamps 
    fid = fopen(fullfile(sessionPathOri,[sessionName,'.res.',spkgs]), 'r');
    Res = fscanf(fid, '%d');
    fclose(fid);
% LOAD spike cluster identities
    fid = fopen(fullfile(sessionPathOri,[sessionName,'.clu.',spkgs]), 'r');
    nClusters = fscanf(fid, '%d', 1);
    Clu = fscanf(fid, '%d');
    fclose(fid);
% LOAD spike waveforms    
    Spk = bload(fullfile(sessionPathOri,[sessionName,'.spk.',spkgs]), [nChannels, inf],0,'short=>single');
    nSpikes = size(Spk, 2)/nSamples;
    Spk = reshape(Spk, [nChannels, nSamples, nSpikes]);


    
% SELECT Artifacs and obvious noise
    gsi = ~(Clu==0|Clu==1);
% REMOVE Artifacs and obvious noise
    Spk = Spk(:,:,gsi);
    Res = Res(gsi);

    disp(['spkTotal: ',num2str(numel(Clu)),'    spkGood: ',num2str(sum(gsi))]);
    
% RESHAPE spike waveform matrix
    nSpikes = size(Spk, 3);
    Spk = reshape(Spk, nChannels,nSpikes*nSamples);

    
% WRITE new spk file
    fid = fopen(fullfile(sessionPath,[sessionName,'.spk.',spkgs]),'w');
    fwrite(fid,Spk,'short');
    fclose(fid);
% WRITE new res file
    fid = fopen(fullfile(sessionPath,[sessionName,'.res.',spkgs]),'w');
    fprintf(fid,'%i\n',Res);
    fclose(fid);
    

    

end
  
system(['ndm_pca ',sessionName,'.xml']);
for spkg = spikeGroupInds,
    features = strrep(num2str(ones(1,numel(Par.SpkGrps(spkg).Channels).*Par.SpkGrps(spkg).nFeatures+1)),' ','');
    system(['nohup KlustaKwik ',...
           sessionName,' ',num2str(spkg),' -MinClusters 50 -MaxClusters 200 ', ...
           '-MaxPossibleClusters 200 -UseFeatures ',features,' &']);
end

% $$$ spkg = 3;
% $$$ fetn = LoadFet(fullfile(sessionPath,   [sessionName,'.fet.',num2str(spkg)]));
% $$$ feto = LoadFet(fullfile(sessionPathOri,[sessionName,'.fet.',num2str(spkg)]));

% $$$ sessionName = 'g10-20130425';
% $$$ sessionPath = fullfile('/storage/gravio/data/processed/nlx/',sessionName);
% $$$ sessionPathOri = fullfile(sessionPath,'ori');
% $$$ fid = fopen(fullfile(sessionPathOri,[sessionName,'.clu.','2']), 'r');
% $$$ nClusters = fscanf(fid, '%d', 1);
% $$$ Clu = fscanf(fid, '%d');
% $$$ nClusters = max(Clu);
% $$$ fclose(fid);
