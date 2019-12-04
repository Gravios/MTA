function [data, timestamps, info] = ecube_load_continuous_data_blocks(filename, varargin)
%load_open_ephys_data_blocks is a function which loads continuous, event, or spike data files into Matlab.
%It is just a modification of the load_open_ephys_data_faster.m (Youngcho Kim) to enable loading of only 
%a chunk of data specified as a range of 1024-samples blocks.
%It allows to convert large sessions into .dat format even with limited RAM.
%NOTE: for now spikes and events are always loaded for the whole file, i.e. BlockRange2Load is ignored in this case.
%
%USAGE:  [data, timestamps, info] = load_open_ephys_data_blocks(filename, <outputFormat>, <BlockRange2Load>)
%
%INPUT:
% filename           is a name of the data file to load.
% <outputFormat>     (optional) If omitted, continuous data is output in double format and is scaled to reflect microvolts. 
%                    If this argument is 'unscaledInt16' and the file contains continuous data, the output data will be in
%                    int16 format and will not be scaled; this data must be manually converted to a floating-point format
%                    and multiplied by info.header.bitVolts to obtain microvolt values. This feature is intended to save memory
%                    for operations involving large amounts of data. Default = [] (double, mkV);
% <BlockRange2Load>  is a 1x2 vector with a range of data blocks (1024 samples each)to be loaded. Default=[] (the whole file).
%
%OUTPUT:
% data:   either an array with continuous samples (in microvolts unless outputFormat is specified, see above),
%         a matrix of spike waveforms (in microvolts),
%         or an array of event channels (integers).
% timestamps:  a vector with timestamps (samples)
% info:        structure with header and other information.
%
%
%EXAMPLE:   [data, timestamps, info] = load_open_ephys_data_blocks(filename, 'unscaledInt16', [25 30])
%
%
%DEPENDENCIES: none
%
% Evgeny Resnik
% version 20.09.2017
%
%
%----------------------------------------------------------------------------------------------------%
%   DISCLAIMER:%
%   Both the Open Ephys data format and this m-file are works in progress.
%   There's no guarantee that they will preserve the integrity of your
%   data. They will both be updated rather frequently, so try to use the
%   most recent version of this file, if possible.
%     ------------------------------------------------------------------%
%     Copyright (C) 2014 Open Ephys
%     ------------------------------------------------------------------
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     <http://www.gnu.org/licenses/>.




% filename = 'TC04-20160609-04-ses4-oephys-100-CH-001.continuous';
% BlockRange2Load = [1 2];



[~,~,filetype] = fileparts(filename);
if ~any(strcmp(filetype,{'.events','.continuous','.spikes'}))
    error('File extension not recognized. Please use a ''.continuous'', ''.spikes'', or ''.events'' file.');
end


%Parse the input parameters
bInt16Out = false;
BlockRange2Load = [];

if nargin > 3
    error('Too many input arguments.');
elseif nargin == 3
    BlockRange2Load = varargin{2};
    if strcmpi(varargin{1}, 'unscaledInt16')
        bInt16Out = true;
    else
        error('Unrecognized output format in "outputFormat".');
    end    
elseif nargin == 2
    if strcmpi(varargin{1}, 'unscaledInt16')
        bInt16Out = true;
    else
        error('Unrecognized output format in "outputFormat".');
    end
end


%Check BlockRange2Load on consistency
if ~isempty(BlockRange2Load)    
    if ~isnumeric(BlockRange2Load) || ~isequal(size(BlockRange2Load(:)) , [2 1])
        error('Input parameter "BlockRange2Load" must be a 1x2 vector.')
    end
    if prod(BlockRange2Load)<=0
        error('Input parameter "BlockRange2Load" must contain only positive values.')
    end
    if BlockRange2Load(1)>=BlockRange2Load(2)
        error('Input parameter "BlockRange2Load" must contain distinct values in ascending order.')
    end
end


fid = fopen(filename);
fseek(fid,0,'eof');
filesize = ftell(fid);

NUM_HEADER_BYTES = 1024;
fseek(fid,0,'bof');
hdr = fread(fid, NUM_HEADER_BYTES, 'char*1');
info = ecube_get_header(hdr);
if isfield(info.header, 'version')
    version = info.header.version;
else
    version = 0.0;
end

switch filetype
    case '.events'
        bStr = {'timestamps' 'sampleNum' 'eventType' 'nodeId' 'eventId' 'data' 'recNum'};
        bTypes = {'int64' 'uint16' 'uint8' 'uint8' 'uint8' 'uint8' 'uint16'};      
        bRepeat = {1 1 1 1 1 1 1};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.2, dblock(7) = [];  end
        if version < 0.1, dblock(1).Types = 'uint64'; end
    case '.continuous'
        SAMPLES_PER_RECORD = 1024;
        bStr = {'ts' 'nsamples' 'recNum' 'data' 'recordMarker'};
        bTypes = {'int64' 'uint16' 'uint16' 'int16' 'uint8'};
        bRepeat = {1 1 1 SAMPLES_PER_RECORD 10};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.2, dblock(3) = []; end
        if version < 0.1, dblock(1).Types = 'uint64'; dblock(2).Types = 'int16'; end
    case '.spikes'
        num_channels = info.header.num_channels;
        num_samples = 40; 
        bStr = {'eventType' 'timestamps' 'timestamps_software' 'source' 'nChannels' 'nSamples' 'sortedId' 'electrodeID' 'channel' 'color' 'pcProj' 'samplingFrequencyHz' 'data' 'gain' 'threshold' 'recordingNumber'};
        bTypes = {'uint8' 'int64' 'int64' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint16' 'uint8' 'float32' 'uint16' 'uint16' 'float32' 'uint16' 'uint16'};
        bRepeat = {1 1 1 1 1 1 1 1 1 3 2 1 num_channels*num_samples num_channels num_channels 1};
        dblock = struct('Repeat',bRepeat,'Types', bTypes,'Str',bStr);
        if version < 0.4,  dblock(7:12) = []; dblock(8).Types = 'uint16'; end
        if version == 0.3, dblock = [dblock(1), struct('Repeat',1,'Types','uint32','Str','ts'), dblock(2:end)]; end
        if version < 0.3, dblock(2) = []; end
        if version < 0.2, dblock(9) = []; end
        if version < 0.1, dblock(2).Types = 'uint64'; end
end

%Structure of a single data block (bytes)
blockBytes = str2double(regexp({dblock.Types},'\d{1,2}$','match', 'once')) ./8 .* cell2mat({dblock.Repeat});

%The total number of blocks in the file
nBlocks = floor((filesize - NUM_HEADER_BYTES)/sum(blockBytes));

%If BlockRange2Load is empty, load the whole file 
if isempty(BlockRange2Load)
    BlockRange2Load = [1 nBlocks];
end
   

if BlockRange2Load(2) > nBlocks
    error(sprintf('Block range exceeds the total number of blocks (%d) in the file',nBlocks))
end


%The number of blocks to load
nBlocks2Load = diff(BlockRange2Load)+1;



switch filetype
    case '.events' %NOT modified by Evgeny Resnik
        timestamps = segRead('timestamps')./info.header.sampleRate;
        info.sampleNum = segRead('sampleNum');
        info.eventType = segRead('eventType');
        info.nodeId = segRead('nodeId');
        info.eventId = segRead('eventId');
        data = segRead('data');
        if version >= 0.2, info.recNum = segRead('recNum'); end              

    
    case '.continuous' %modified by Evgeny Resnik
        if nargout>1
            %load timestamps of data blocks
            info.ts = segRead_blocks('ts', BlockRange2Load, fid, dblock, NUM_HEADER_BYTES, blockBytes);
            
%             %Load timestamps of data blocks
%             %index of the data type 
%             segNum = find(strcmp({dblock.Str},'ts'));            
%             %set file position indicator
%             start_pos = NUM_HEADER_BYTES + (BlockRange2Load(1)-1)*sum(blockBytes) + sum(blockBytes(1:segNum-1));
%             fseek(fid, start_pos, 'bof');            
%             %load timestamps of blocks (!) from the specified range of blocks (int64 into double)
%             nSamples2Load = nBlocks2Load*dblock(segNum).Repeat;
%             Precision     = sprintf('%d*%s', dblock(segNum).Repeat,dblock(segNum).Types);
%             Skip          = sum(blockBytes) - blockBytes(segNum); 
%             info.ts2 = fread(fid, nSamples2Load, Precision, Skip, 'l'); 
%             clear segNum start_pos nSamples2Load Precision Skip            
        end
        
        %load nsamples info for the blocks
        info.nsamples = segRead_blocks('nsamples', BlockRange2Load, fid, dblock, NUM_HEADER_BYTES, blockBytes);

        
%         %index of the data type
%         segNum = find(strcmp({dblock.Str},'nsamples'));
%         %set file position indicator
%         start_pos = NUM_HEADER_BYTES + (BlockRange2Load(1)-1)*sum(blockBytes) + sum(blockBytes(1:segNum-1));
%         fseek(fid, start_pos, 'bof');
%         %load info from the specified range of blocks (int64 into double)
%         nSamples2Load = nBlocks2Load*dblock(segNum).Repeat;
%         Precision     = sprintf('%d*%s', dblock(segNum).Repeat,dblock(segNum).Types);
%         Skip          = sum(blockBytes) - blockBytes(segNum);
%         info.nsamples = fread(fid, nSamples2Load, Precision, Skip, 'l');
%         clear segNum start_pos nSamples2Load Precision Skip
        
        
        if ~all(info.nsamples == SAMPLES_PER_RECORD)&& version >= 0.1, error('Found corrupted record'); end
        if version >= 0.2, 
            info.recNum = segRead_blocks('recNum', BlockRange2Load, fid, dblock, NUM_HEADER_BYTES, blockBytes);
        end
        
        %load data blocks
        data = segRead_blocks('data', BlockRange2Load, fid, dblock, NUM_HEADER_BYTES, blockBytes);

        
%         %index of the data type 
%         segNum = find(strcmp({dblock.Str},'data'));        
%         %set file position indicator
%         %block structure: 12 bytes - not used, 1024 bytes - data, 10-bytes - not used       
%         start_pos = NUM_HEADER_BYTES + (BlockRange2Load(1)-1)*sum(blockBytes) + sum(blockBytes(1:segNum-1));
%         fseek(fid, start_pos, 'bof');
%         %load data from the specified range of blocks (keep as int16)
%         nSamples2Load = nBlocks2Load*dblock(segNum).Repeat;
%         Precision = [sprintf('%d*%s', dblock(segNum).Repeat,dblock(segNum).Types) '=>int16'];
%         Skip = sum(blockBytes) - blockBytes(segNum);    
%         data2 = fread(fid, nSamples2Load, Precision, Skip, 'b');
%         clear segNum start_pos nSamples2Load Precision Skip      
    
        %Create timestamp vector (samples) for all samples if requested
        if nargout>1
            timestamps = nan(size(data));
            current_sample = 0;
            for record = 1:length(info.ts)
                timestamps(current_sample+1:current_sample+info.nsamples(record)) = info.ts(record):info.ts(record)+info.nsamples(record)-1;
                current_sample = current_sample + info.nsamples(record);
            end
        end
        

        %optional convertion to mkV (from int16 to double)
        if ~bInt16Out
            data = double(data) .* info.header.bitVolts;
        end
        

        
        
        
    case '.spikes' %NOT modified by Evgeny Resnik
        timestamps = segRead('timestamps')./info.header.sampleRate;
        info.source = segRead('source');
        info.samplenum = segRead('nSamples');
        info.gain = permute(reshape(segRead('gain'), num_channels, numIdx), [2 1]);
        info.thresh = permute(reshape(segRead('threshold'), num_channels, numIdx), [2 1]);
        if version >= 0.4, info.sortedId = segRead('sortedId'); end
        if version >= 0.2, info.recNum = segRead('recordingNumber'); end
        data = permute(reshape(segRead('data'), num_samples, num_channels, numIdx), [3 1 2]);
        data = (data-32768)./ permute(repmat(info.gain/1000,[1 1 num_samples]), [1 3 2]);
end

fclose(fid);




end %end of function


%=============================== Supplementary functions =========================================%


function seg = segRead_int16(segName, mf)
%This function is specifically for reading continuous data.
%  It keeps the data in int16 precision, which can drastically decrease
%  memory consumption
if nargin == 1, mf = 'l'; end
segNum = find(strcmp({dblock.Str},segName));
fseek(fid, sum(blockBytes(1:segNum-1))+NUM_HEADER_BYTES, 'bof');
seg = fread(fid, numIdx*dblock(segNum).Repeat, [sprintf('%d*%s', ... 
    dblock(segNum).Repeat,dblock(segNum).Types) '=>int16'], sum(blockBytes) - blockBytes(segNum), mf);
end


function seg = segRead(segName, mf)
    if nargin == 1, mf = 'l'; end
    segNum = find(strcmp({dblock.Str},segName));
    fseek(fid, sum(blockBytes(1:segNum-1))+NUM_HEADER_BYTES, 'bof'); 
    seg = fread(fid, numIdx*dblock(segNum).Repeat, sprintf('%d*%s', ...
        dblock(segNum).Repeat,dblock(segNum).Types), sum(blockBytes) - blockBytes(segNum), mf);    
end



function seg = segRead_blocks(segName, BlockRange2Load, fid, dblock, NUM_HEADER_BYTES, blockBytes)
%Supplementary function for load_open_ephys_data_blocks.m
% Evgeny Resnik
% version 20.09.2017

%The number of blocks to load
nBlocks2Load = diff(BlockRange2Load)+1;

%index of the data type
segNum = find(strcmp({dblock.Str}, segName));

nSamples2Load = nBlocks2Load*dblock(segNum).Repeat;
Skip = sum(blockBytes) - blockBytes(segNum);

%set file position indicator
start_pos = NUM_HEADER_BYTES + (BlockRange2Load(1)-1)*sum(blockBytes) + sum(blockBytes(1:segNum-1));
fseek(fid, start_pos, 'bof');

switch segName
    case {'ts', 'nsamples', 'recNum'}
        Precision     = sprintf('%d*%s', dblock(segNum).Repeat,dblock(segNum).Types);
        ml = 'l';

    case 'data'
        Precision = [sprintf('%d*%s', dblock(segNum).Repeat,dblock(segNum).Types) '=>int16'];
        ml = 'b';       
        
    otherwise
        error('Unknown segName value!')
end

seg = fread(fid, nSamples2Load, Precision, Skip, ml);

end %segRead_blocks


