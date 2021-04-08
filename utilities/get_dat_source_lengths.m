function srsLength = get_dat_source_lengths(meta)
% function srsLength = get_dat_source_lengths(meta)
%
% Return the length in system samples of all parts which were concatenated into the dat file
%
% ARGIN:
%     meta - struct: contains paths to the files
%            required fields: 
%                meta.path.processed.ephys
%                meta.sessionName
%

Par = LoadPar(fullfile(meta.path.processed.ephys, [meta.sessionName '.xml']));            
if exist(fullfile(meta.path.processed.ephys,[meta.sessionName,'.srs']))
    srsLength = parse_dat_source_lengths(meta);
else
    f = dir(fullfile(meta.path.processed.ephys,[meta.sessionName,'.dat']));
    srsLength = {f.bytes/Par.nChannels/(Par.nBits/8)};
end
