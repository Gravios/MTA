function TTLValue = convert_srs_to_TTLValue(meta)
% function TTLValue = convert_srs_to_TTLValue(meta)
% 
% generate mapping of recording to mocap
% 
% ARGIN:
%     meta - struct: 
%         required fields:
%             meta.sessionName
%             meta.csv
%             meta.path.processed.ephys
%


% GET record indices where meta.csv is occupied
csvIndex = find(~cellfun(@isempty,meta.csv));
% GET the length of each record in primary acquisiton system sampling rate
recordLengths = get_dat_source_lengths(meta);
% GENERATE take list and prime subject list which matches the meta.csv list
takeList = repmat({''},size(meta.csv));
primeSubject = repmat({''},size(meta.csv));
for t = 1:numel(csvIndex),
    takeList{csvIndex(t)} = [meta.sessionName ,num2str(t,'.take_%04.f.mat')];
    primeSubject{csvIndex(t)} = meta.primarySubject;
end

% CONCATENATE into sudo TTLValue list
TTLValue = cf(@(r,t,p) {r,t,p}, recordLengths, takeList, primeSubject);
