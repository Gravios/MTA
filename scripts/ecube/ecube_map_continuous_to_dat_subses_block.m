function ecube_map_continuous_to_dat_subses_block( filebase, subSesFileNames, dat, chanInfo, varargin)
% function map_oephys_to_dat_subses_block(subSesFileNames, filebase, dat, acqSystem, chanInfo, datFileSizeLimit)
% oephys2dat_subses_blocks is a modification of oephys2dat_subses.m, which
% converts .continuous files from "subSesFileNames" with wideband signal
% recorded with an open ephys acquisition system to .dat data format, but
% does it in small blocks. It allows to deal with large sessions.
% This function is a supplementary function for oephys2dat and used to convert files within 
% individual subsessions (sleep_box, maze, arena and so on) of the experiment.
% Use oephys2dat.m to convert the whole experiment session.
%
% NOTE : although CH, ADC and AUX channels have different header.bitVols values, for simplicity only the CH-value (0.195) is applied to all
%      channels to convert to microvolts.Take it into account if you need real amplitudes on ADC/AUX channels.
%
% USAGE : oephys2dat_subses_blocks(subSesFileNames, dat, <datFileSizeLimit>)
%
% INPUT :
% subSesFileNames       is a cell vector with a list of .continuous files to be converted.
%                Files must be generated by the same open ephys processor and belong to the same subsession. 
% dat        is a name of the output .dat file the converted data will be saved in.
%                Conventional naming scheme: filebase-subSesId-subSesDescription.dat 
% acqSystem      is a name of the aqcuisition system ('oephys' - open ephys; 'ecube' - eCube; 'nlx' - Neuralynx).
% <datFileSizeLimit>  is the maximum size of simultaneously loaded to RAM data. Default = 40 (empirically tested on the labs erver).
%
% OUTPUT :
% FileOut.dat      is a binary .dat file with the wideband signals (mkV).
% FileOut.dat.info is an ASCII file with info about FileOut.dat [nChannels nSamples SamplingRate_Hz].
% FileOut.dat.sts  is an ASCII file with original timestamps for samples in FileOut.dat (.dat samples)
%
% EXAMPLE : oephys2dat_subses_blocks(subSesFileNames, 'ER50-20170101-01-hc.dat', 'oephys')
%
% DEPENDENCIES :
% labbox: DefaultArgs, memorylinux, start_parallel
% oephys: load_open_ephys_header
% FMAToolbox: SaveBinary_FMAToolbox
%
% Evgeny Resnik
% version 12.04.2016
% version 07.06.2018: can deal with eCube data files
%
% Gerrit
% version 18.02.2019: some changes to deal with recording
% gaps and channel groups of different start/stop times, data loss or cummulative channel misalignment should be gone now, also shortened
% code


% DEFARGS ------------------------------------------------------------------------------------------
[ acqSystem, datFileSizeLimit ] = DefaultArgs(varargin,{ 'ecube', 40 });
source = [filebase,'-',acqSystem];
nFiles = length(subSesFileNames);
%---------------------------------------------------------------------------------------------------

% CHECK that all the files exist and that file names and internal header info match
for k=1:nFiles    
    if ~exist(subSesFileNames{k}, 'file')
        error(sprintf('map_oephys_to_dat_subses_block: File %s not found', subSesFileNames{k}))
    end
    
    header(k) = ecube_load_header(subSesFileNames{k});    
    %Deal with a files from eCube with additional suffixes "HS_" in file name, get rid of it first
    header(k).channel = strrep(header(k).channel,'HS_','');

% REMOVE LATER
    chanInfo(k).chanId = regexprep(chanInfo(k).chanId,'(^0*)(\d+)','$2');
% LATER REMOVE

    assert(strcmp(header(k).channel,[chanInfo(k).chanType,chanInfo(k).chanId]),...
           'ecube_map_continuous_to_dat_subses_block:HeaderChannelIdMismatch');
end%for k


%%%<<< CALCULATE parameters for a loop across data chunks

%Compute size of files in the list (bytes)
for k=1:length(subSesFileNames);
    fid = fopen(subSesFileNames{k});
    fseek(fid,0,'eof');
    fileSizeBytes(k) = ftell(fid);
    fclose(fid);
end

% SET structure of a block of .continuous data
blockBytes =  [8 2 2 2048 10]; 

% COMPUTE size of a data block (bytes)
blockSamples = 1024;
% SET Pre-defined size of a .continuous file header (bytes)
NUM_HEADER_BYTES = 1024;

% COMPUTE the total number of blocks in the file
nBlocks = (fileSizeBytes - NUM_HEADER_BYTES)/sum(blockBytes);
assert(all(floor(nBlocks)==nBlocks),[mfilename,':IncompleteBlock'])

% SET LENGTH of a data chunk (in blocks) that can be loaded simulateously across channels
chunkBlock = 512;
chunkSamples = chunkBlock*blockSamples;

% COMPUTE number of full chunks in the file
nChunks = ceil(max(nBlocks/chunkBlock));
firstTimeStamp = nan([size(subSesFileNames)]);
finalTimeStamp = nan([size(subSesFileNames)]);
for k=1:nFiles
    [~, t, ~] = ecube_load_continuous_data_blocks(subSesFileNames{k}, ...
                                            'unscaledInt16', ...
                                            [1,2]);
    firstTimeStamp(k)=min(t);
    [~, t, ~] = ecube_load_continuous_data_blocks(subSesFileNames{k}, ...
                                            'unscaledInt16', ...
                                            [nBlocks(k)-1 nBlocks(k)]);
    finalTimeStamp(k)=max(t);
end
clear('t');

startTimeStamp = min(firstTimeStamp);

%%%>>>


%%%<<< Load from all channels and save to a .dat file
% Does not handle gaps in the data

% LOAD / CHECK / CORRECT / SAVE all chunks

% SET startblock
blockRange = repmat([1   chunkBlock],nFiles,1)';    
if blockRange(2,1) > min(nBlocks),  blockRange(2,:) = min(nBlocks);  end
remnantData       = repmat({int16([])},[1,nFiles]);
remnantTimeStamps = repmat({[]},[1,nFiles]);
finished=false;
n=0;    
finalDat = [];

while ~finished;
    n=n+1;
    
    if any(blockRange(2,:)>nBlocks),
        blockRange(2,blockRange(2,:)>nBlocks) = nBlocks;
    end
    
% STOP when all data blocks from all files have been read      
    if any(blockRange(2,:)==nBlocks), finished = true; end
    
% LOADING section: exactly as in oephys2dat_subses.m, but for a chunk
    %Initialize the variable for accumulating the data
    %Note: a cell array is used here to be able to run parfor and to accomodate cases when channels have different length
    rawData   = cell([1,nFiles]);
    rawTimeStamps = cell([1,nFiles]); 
    
    for k=1:nFiles
% LOAD chunk
        [rawData{k}, rawTimeStamps{k}] = ...
            ecube_load_continuous_data_blocks(subSesFileNames{k},'unscaledInt16',blockRange(:,k));
% STORE data and timestamps
        if diff(blockRange(:,k))<=1;
            rawData{k}=0;
            t=0;
        end
    end%for k
    assert(all(cellfun(@(x) all(diff(x)==1), rawTimeStamps)),[mfilename,':TimeStampGap']);

    rawData       = cf(@(x,y) cat(1,x,y), remnantData,rawData);
    rawTimeStamps = cf(@(x,y) cat(1,x,y), remnantTimeStamps,rawTimeStamps);    

    
% FILL data into the right timeslots

    if n==1,
        fid = fopen(dat,'a');        
        timeStampShift = startTimeStamp;
    else
        timeStampShift = chunkLastTimeStamp+1;
    end

    chunkStartTimeStamp  = min(cellfun(@(x) x(1), rawTimeStamps));
    chunkLastTimeStamp  = min(cellfun(@(x) x(end), rawTimeStamps));
    chunkRange = chunkLastTimeStamp-chunkStartTimeStamp+1;
    
    assert( timeStampShift == chunkStartTimeStamp ,[mfilename,':GapBetweenBlocks']);
    
% INITIATE new raw vector        
    processedData = repmat({zeros(1,chunkRange)},[nFiles,1]); 
    for k=1:nFiles
        shiftedTimeStamps{k} = rawTimeStamps{k}-timeStampShift+1;
        remnantTimeInds{k}   = shiftedTimeStamps{k} >  chunkRange;
        remnantTimeStamps{k} = rawTimeStamps{k}(remnantTimeInds{k});
        shiftedTimeStamps{k} = shiftedTimeStamps{k}(shiftedTimeStamps{k} <= chunkRange);
    end
        
    for k=1:nFiles
        processedData{k}(shiftedTimeStamps{k}) = rawData{k}(shiftedTimeStamps{k});
    end
    
    for k=1:nFiles
        if ~isempty(remnantTimeStamps{k}),
            remnantData{k} = rawData{k}(remnantTimeInds{k});
        end
    end
    
    assert(any(chunkLastTimeStamp+1 == min(cellfun(@(x) x(1), ...
               remnantTimeStamps(~cellfun(@isempty,remnantTimeStamps))))),...
           [mfilename,':GapBetweenBlocks']);

% SAVE / APPEND data (mkV) to .dat file
    fprintf('Saving data chunk-%d of %d (1024-samples, %d channels) to %s ... ', ...
            n, nChunks, nFiles, dat);
    fwrite(fid, cell2mat(processedData),'int16');
    fprintf('DONE\n')
    clear('rawData', 'rawTimeStamps', 'processedData', 'processedTimeStamps');

    
% END of Loading section ----------------------------------------------------%

% $$$    finalDat = cat(1,finalDat,cell2mat(processedData)');
   blockRange = blockRange + chunkBlock;      
   
end% loop across chunks

fclose(fid);

%%%>>>

% CREATE an info file: [nChan sampleRate nSamples StartTime EndTime]

%% NOTE : if the data contained time gaps, a new time vector without gaps is created and 
%% StartTime/ EndTime are taken from this new time vector. 
dat_info = [dat '.info'];
fprintf('Saving [nChan sampleRate nSamples startTime endTime] into a file %s ...', dat_info)
info = [ numel(header),                       ... nChan
         header(1).sampleRate,                ... sampleRate
         nBlocks*1024,                        ... nSamples
         min(firstTimeStamp),                 ... start time
         min(firstTimeStamp)+nBlocks*1024-1   ... end time    
       ]; 
dlmwrite(dat_info, info, 'delimiter','\n','precision','%.0f');
fprintf('DONE\n')



return









% DIAGNOSTIC figures
% $$$ figure;plot(fileSizeBytes)
% $$$ xlabel('files')
% $$$ ylabel('File Size bytes')
% $$$ title([cdir ' - FileSize_bytes for each file - '])
