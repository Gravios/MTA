function repair_ncs_timestamp_order(varargin)

% DEFARGS ----------------------------------------------------------------------------------%    
defargs = struct (...
    'ncsPath',           '/storage/gravio/data/raw/nlx/jg05-20120312-nlx/',...
    'ncsFilebase',       'jg05-20120312-01-cof-nlx-'                      ,...
    'nChannels',         96                                               ,...
    'NCS_EXT',           '.ncs'                                           ,...      
    'HEADER_SIZE',       16*1024                                          ,... % bytes
    'SAMPLE_RATE',       32556                                            ,... % Hertz
    'DBLOCK_HEADER_SIZE',8+4+4+4                                          ,... % bytes - uint64 + uint32(3)
    'DATA_SIZE',         2);                                                   % bytes - uint16

[ncsPath, ncsFilebase, nChannels, NCS_EXT, HEADER_SIZE, SAMPLE_RATE,...
DBLOCK_HEADER_SIZE, DATA_SIZE] = DefaultArgs(varargin,defargs,'--struct');
% END DEFARGS ----------------------------------------------------------------------------------%    

keyboard

% MAIN ----------------------------------------------------------------------------------%

cpwd = pwd;  % Save path of current directory
cd(ncsPath); % Move workspace to ncsPath 

% COLLECT information from all records in each file
blockTimeStamp  = [];
blockChan       = [];
blockSampleRate = [];
blockSampleCount= [];

BlockChansAgree = false([1,nChannels]);  
BlockSampleRatesAgree = false([1,nChannels]);

for chan = 1:nChannels,

    godNcsFid = fopen(fullfile(ncsPath,[ncsFilebase,sprintf('%03d',chan),NCS_EXT]),'a+');

    fseek(godNcsFid, 0, 'eof');
    godEof = ftell(godNcsFid);

    fseek(godNcsFid, 0, 'bof');    
    HeaderG = fread(godNcsFid, HEADER_SIZE, '*char');

    b = 1;
    while ftell(godNcsFid) < godEof,
        blockTimeStamp(b,chan)  = fread(godNcsFid, 1, '*uint64'); % block timestamp
        blockChan(b,chan)       = fread(godNcsFid, 1, '*uint32'); % block channel
        blockSampleRate(b,chan) = fread(godNcsFid, 1, '*uint32'); % block SAMPLE_RATE
        blockSampleCount(b,chan)= fread(godNcsFid, 1, '*uint32'); % block samples
        fseek(godNcsFid, blockSampleCount(b,chan)*2, 'cof');      % block Data
        b = b+1;
    end
    b = b-1;
    

    if all(blockSampleRate(1:b,chan)==SAMPLE_RATE), BlockSampleRatesAgree(chan) = true; end
    if all(blockChan(1:b,chan)==chan-1), BlockChansAgree(chan) = true; end
    
    fclose(godNcsFid);
end

tsFinal = zeros([1,nChannels]);
tsStart = zeros([1,nChannels]);
bcTs = {};
for chan = 1:nChannels,
    bcTs{chan} = blockTimeStamp(blockTimeStamp(:,chan)~=0,chan);
    tsStart(chan) = bcTs{chan}(1);
    tsFinal(chan) = bcTs{chan}(end);        
end

FirstTimesStampMatch = false; 
if all(tsStart == tsStart(1)), 
    FirstTimesStampMatch = true; 
end

FinalTimesStampMatch = false; 
if all(tsFinal == tsFinal(1)), 
    FinalTimesStampMatch = true; 
end

degenRecordCount = cell2mat(cellfun(@(x) numel(x)-numel(unique(x)),bcTs,'UniformOutput',false));
NoDegenerateTimestamps = false;    
if ~any(degenRecordCount)
    NoDegenerateTimestamps = true;
end
    

% Assert some assumptions are actually true.
assert(FirstTimesStampMatch,'MTA:utilities:repair_ncs_timestamp_order:FirstTimesStampDoNotMatch');
assert(FinalTimesStampMatch,'MTA:utilities:repair_ncs_timestamp_order:FinalTimesStampDoNotMatch');


% Sort timestamps and get new write order
[sbcTs,sbcTsInds] = cellfun(@sort,bcTs,'UniformOutput',false); 

degenRecordTimestamps = [];
if ~NoDegenerateTimestamps,
    for chan = find(degenRecordCount)
        degenRecordTimestamps = [degenRecordTimestamps;sbcTs{chan}(find(diff(sbcTs{chan})==0))];
    end
end

ats = [];
for chan = 1:nChannels, ats = [ats;sbcTs{chan}]; end
ats = unique(ats);

% remove untrustworthy records with degenereate timestamps
ats(ismember(ats,degenRecordTimestamps))=[];


% @ 39min 18 sec a gap occurs which causes a misalignment

sampleRate = uint32(SAMPLE_RATE);
nSamples =   uint32(512);

my_ats = ats;
bad_ats = find(diff(my_ats)<1e4);
my_bad_ats = my_ats(bad_ats);
my_ats(bad_ats) = [];
%my_ats([bad_ats(1:2:end-1)-1;bad_ats(end)+1]) = [];    
%cind = 2500000:numel(my_ats);
%bts = my_ats(cind);
bts = my_ats;
for chan = 1:nChannels,
    tic
    % Open original NCS file.
    godNcsFid = fopen(fullfile(ncsPath,[ncsFilebase,sprintf('%03d',chan),NCS_EXT]),'r');
    % Open repair NCS file.
    repNcsFid = fopen(fullfile(ncsPath,[ncsFilebase,sprintf('reorder-%03d',chan),NCS_EXT]),'w+'); 

    HeaderG = fread(godNcsFid, HEADER_SIZE/2, '*int16');    % Get Header from original file
    fwrite(repNcsFid, HeaderG, 'int16');                    % Write Header to repair file

    chanId = uint32(chan-1);

    tind = sbcTs{chan}>=bts(1)&bts(end)>=sbcTs{chan};
    tts = sbcTs{chan}(tind);
    bad_bts = ismember(tts,my_bad_ats);
    tts(bad_bts) = [];
    chunk = sbcTsInds{chan}(tind);
    chunk(bad_bts) = [];
    mits = fliplr(bts(~ismember(bts,tts))');
    
    for b = 1:numel(chunk), % FE data block, write in correct order 
        
        fseek(godNcsFid,...                               % Skip the ...
              HEADER_SIZE+...                             %     file header,
              DBLOCK_HEADER_SIZE.*chunk(b)+... %     all preceding data block headers,
              DATA_SIZE.*512.*(chunk(b)-1),... 
              'bof');                                     % starting frnom the begining of the file
        
        if ~isempty(mits),
            while tts(b)>mits(end), % write filler record
                write_dummy_record(repNcsFid,mits(end),chanId,sampleRate,nSamples)
                mits(end) = [];
                if isempty(mits), break; end
            end
        end

        % write data block header
        fwrite(repNcsFid, uint64(tts(b)), 'uint64');     % block timestamp
        fwrite(repNcsFid, chanId,         'uint32');     % block channel
        fwrite(repNcsFid, sampleRate,     'uint32');     % block SAMPLE_RATE
        fwrite(repNcsFid, nSamples,       'uint32');     % block samples            
                                                    % read data from original file at sorted index b write data to repair file
        fwrite(repNcsFid,fread(godNcsFid,512,'*uint16'),'uint16');

    end

    while ~isempty(mits), % write filler record
        write_dummy_record(repNcsFid,mits(end),chanId,sampleRate,nSamples)
        mits(end) = [];
    end
    
    fclose(godNcsFid); % close original file
    fclose(repNcsFid); % close repair file
    toc
end

cd(cpwd); % Return to original directory

% END MAIN ----------------------------------------------------------------------------------%



