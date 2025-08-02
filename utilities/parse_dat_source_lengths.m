function srsLength = parse_dat_source_lengths(meta)
% function srsLength = parse_dat_source_lengths(meta)
% 
% parse frome .srs file the length of each 
%
% each line is the source file and data point count separated by a space character
%
fPosition = fullfile(meta.path.processed.ephys,[meta.sessionName,'.srs']);
fid     = fopen(fPosition,'r');

srsLength = {};
l = 1;
srsline = fgetl(fid);
while srsline ~= -1
    srsline = strsplit(srsline,' ');
    srsLength{l} = str2num(srsline{2});
    srsline = fgetl(fid);
    l = l + 1;
end
fclose(fid);

