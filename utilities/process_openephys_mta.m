function process_openephys_mta(ephysSource,ephysDestination,subjectName,sessionName,recordingTypes)
% REQVARS for processing raw openephys directories 


% TESTING vars


% $$$ subjectName = 'RS0317';
% $$$ sessionName = 'RS0317-20170517';
% $$$ recordingTypes = {'lfp','fbr'};
% $$$ ephysSource      = fullfile('/storage/gravio/data/raw/openephys/',subjectName,[sessionName,'-openephys']');
% $$$ ephysDestination = '/storage/gravio/data/processed/nlx/';

sessionEphysPath = fullfile(ephysDestination,subjectName,sessionName);

try,mkdir(sessionEphysPath);end
processedOpenephysFilebase = fullfile(ephysDestination,subjectName,sessionName,sessionName);

subSessionDir = dir(ephysSource);
subSessionDir(1:2) = [];
subSessionDir(~[subSessionDir.isdir]) = [];
copyfile(fullfile(ephysSource,subjectName,[sessionName,'-openephys'],...
                  subSessionDir(1).name,[subSessionDir(1).name,'.xml']),...
         [ephysDestination,'.xml']);

Par = LoadPar([ephysDestination,'.xml']);

for rt = recordingTypes
    system(['touch ',processedOpenephysFilebase,'.' rt{1}])
end

system(['ls -lh ',sessionEphysPath])

% Concatenate lfp files 
s = 1;
for subSession = subSessionDir',
    subsessionFilebase = fullfile(ephysSource,subjectName,[sessionName,'-openephys'],subSession.name,subSession.name);
    system(['cat ',subsessionFilebase,'.lfp >> ',processedOpenephysFilebase,'.lfp']);
    lfpSampleCount(s) = numel(LoadBinary([subsessionFilebase,'.lfp'],1,Par.nChannels,4));
    fbr(s) = load([subsessionFilebase,'.fbr'],'-mat');   
    fbr(s).photo_time = [lfpSampleCount(s)];
    s = s+1;
end            

save([processedOpenephysFilebase,'.fbr.mat'],'fbr');