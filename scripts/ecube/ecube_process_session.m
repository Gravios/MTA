function ecube_process_session(filebase, xml, varargin)
% function  process_raw_ecube(filebase, xml, <channelMap>, <processorList>)
%
% <UserNAME>    : regexp : [a-z][A-Z]+
% <UserID>      : regexp : [a-z][A-Z]{1,4}
% <AnimalID>    : regexp : [a-z][A-Z]{1,4}\d{1,4}
%
% TEMPLATES :
%     metadata:   oephys/template_files/
%     parameters: oephys/template_files/
%
% PREPARE data for processing:
%    
%    PLACE Place all subsessions in subject directory in /storage/<User>/data/raw/ecube/<AnimalID>/.
%    CREATE an ASCII *.note file with experiment's metadata. %            
%    PLACE *.note in processing directory. 
%    RENAME *.note template copy as AnimalID-YYYYMMDD.note 
%    UPDATE *.note with experiment's metadata
%    CREATE XML parameter file for data processing. (AKA "parameter file")   
%    UPDATE the following fields within parameter file
%
%        General information:
%            Description: description1,description2,… a list of single-word descriptions 
%                         for subsessions of an experiment (e.g. sleep,maze,sleep or ca1, ca3).
%
%        Acquisition System:
%            Number of channels: total number of channels recorded, including AUX and ADC channels.
%            Sampling rate: sampling rate of the acquisition system (default = 30000).
%            Initial offset: 0
%            Amplification: 1 (for standard intant headstages without on-board amplification).
%
%        Local Field Potentials:
%            Sampling rate: sampling rate of downsampled *.lfp file. default=1250. 
%
%        Files:
%        Extension: lfp
%        Sampling rate: 1250
%     
%        Channel mapping": a list of open ephys channel indicies in a proper anatomical order. 
%                          NOT USED AT THE MOMENT
%
%        Anatomical Groups:
%            Groups: sets of channel indicies grouped by electrode/probe site configuration.
%                    Do not forget to put all accelerometer (AUX) channels into a separate anatomical group.
%
%        Spike Groups:
%            Groups": groups of open ephys channel indicies grouped in accordance with appearance of
%                     spikes from the same neurons. 
%                     It will not be needed for KiloSort. 
%                     BUT NEEDED FOR OTHER ANALYSIS
% 
%    CREATE an ASCII .chanmap file 
%        a list of open ephys channel indicies in a proper anatomical order [ 1 9,4,3,...]
%            This file will be used by oephys2dat.m to reorder channels in the *.dat file. 
%            The mapping file must contain ONLY indices of LFP channels (CH, not AUX, ADC and so on). 
%            {oephys2dat.m} sorts all the LFP channels and put them first in the created .dat file, 
%            {oephys2dat.m} then sorts channels within each other (AUX, ADC and so on) 
%                           group and append them to the .dat file. 
%
%    IMPORTANT : Preliminary testing has revealed
%        {oephys2dat.m} may crash without showing the "Out of memory" error messages when 
%           running on local PCs. In this case run MATLAB on the lab server.
%
%
% process_raw_ecube is a function which preprocess a session recorded with an open ephys aqcuisition system:
%     LINK raw data to processing directory:  (bash)       {ecube_link_files}
%     RENAME raw data to ndm_naming schema:   (bash)       {ecube_rename_files}
%     CONVERT .continous files to .dat file:  (matlab)     {oephys2dat.m}
%     MERGE .dat files of subsessions:        (?)
%     DOWNSAMPLE .dat file to .lfp (1250Hz):  (c)          {ndm_lfp}
%     IF spike groups are defined in an .xml file)
%         FILTER .dat highpass to .fil:           (c)      {ndm_hipass}
%         DETECT spikes in .fil:                  (c)      {ndm_extractspikes}
%         EXTRACT spikes from .dat file:          (c)      {ndm_extractspikes}
%         COMPUTE their PC features (only  ndm_pca
%        
% Dependencies :
%    NDM toolbox (L. Hazan, M. Zugaro, http://neurosuite.sourceforge.net).
%
%    NOTE: This function can deal with the case when .continuous files from some of the channels are not present. 
%      e.g. when one records with a 32-channels headstage and then delete some channel files.
%      {oephys2dat.m} detects all present .continuous files
%      {oephys2dat.m} groups them by their type (CH, AUX, ADC) 
%      {oephys2dat.m} sorts the channels within each group in ascending order using their channel indices.
% 
% AFTER process_raw_ecube:
%    CHECK filebase.dat and filebase.lfp files in Neuroscope:
%       -mark permanently bad channels as 'skipped'
%       -identify a channel with strong theta oscillation to be used for detecting 
%        brain states (RUN, SWS, REM), i.e. CA1pyr channel.
% 
%    CREATE an ASCII file with additional parameters filebase.mypar
%        NOTE: In case when a recording session consists of a single behavioral episode and the 
%              whole session can be assigned to a particular type of behavior (e.g. RUN or SLEEP),
%              use the following format filebase.mypar:
%            "SLEEP   filebase" instead of standard "SLEEP  hc" (without quotation marks).
% 
% USAGE : process_raw_ecube(filebase, xml, <channelMap>, <processorList>)
%
% INPUT :
%     filebase,   string: name of the session to be processed.
%
%     xml,        string: XML file with parameters necessary for processing
%
%     channelMap, string: Name of ASCII file, which contains list of channel indicies in 
%                 anatomical order (raw- or columnwise).
%                 IF not present, channel mapping is skipped.
%
%     processorList, numericArray: ID of the open ephys recording processor from which the data
%                    files are merged. Default=[];
%                    (e.g. a file 100_CH1.continuous was recorded by the processor 100).
%
% OUTPUT :
%     The script creates a new directory AnimalID-YYYMMDD next to the directory with the raw data 
%     fiels AnimalID-YYYMMDD-oephys and saves there all the converted files: AnimalID-YYYMMDD.dat, 
%     AnimalID-YYYMMDD.lfp, AnimalID-YYYMMDD.cat.evt and additional ASCII .info files with info 
%     about processed subsessions). All processing steps are logged into an ASCII file 
%     AnimalID-YYYMMDD-oephys.ProcessOephys.log.
%
% EXAMPLE :   
%     ProcessEcube('TestEcube-20180619-ecube', '1ses-160chan.xml', [])
%
% DEPENDENCIES :
%     labbox: DefaultArgs, LoadXml 
%     oephys: CorrectNames, oephys2dat, CleanProcessedSession
%     NDM toolbox (not matlab): ndm_lfp, ndm_hipass, ndm_extractspikes, ndm_pca
% 
% Evgeny Resnik
% version 10.09.2017
% Justin Graboski
% version 29.07.2019 refactoring
% 
% TO DO:
%     - fix incorrect renaming of non-default files in the session folder (e.g. avi).
%     - oephys2dat_subses_blocks - replace cell2mat(raw_t) with a loop as in oephys2dat_subses. 
%     - oephys2dat_subses_blocks - make proper calculation of start/end times for time gaps if 
%                              present in the data. 
%     - convert open ephys event file
%     - processing of bonsai video files
%

% DEFARGS ------------------------------------------------------------------------------------------

if nargin<1
    error(['USAGE:  ProcessEcube(Session, xml, <channelMap>, <processorList>)'])
end

[acqSystem, channelMap, processorList, overwrite] = DefaultArgs(varargin,{ 'ecube', [], [], false });

%---------------------------------------------------------------------------------------------------


%%%<<< TESTING VARS
addpath('/storage/gerrit/code/matlab/Wraps/Convert/openephys/')
addpath('/storage/gerrit/code/matlab/Wraps/Convert/Ecube/')

% $$$ filebase      = 'IF13-20190409c';
% $$$ xml           = 'IF13.xml';
% $$$ chanmap    = 'IF13.chanmap';
% $$$ processorList = [100,103]; 
% $$$ acqSystem     = 'ecube';
% $$$ mfilename     = 'process_raw_ecube';

filebase      = 'IF13-20190409b';
xml           = 'IF13.xml';
chanmap       = 'IF13.chanmap';
processorList = [100,103]; 
acqSystem     = 'ecube';
mfilename     = 'process_raw_ecube';

%%%>>>

%%%<<< Preprocessing
source = [filebase,'-',acqSystem];
fprintf(['%s - %s\n'], source, mfilename)
subject = regexprep(filebase,'(\w+)-(\w+)','$1');
create_directory(filebase);

% START logging the command window messages
fileOutLog = sprintf('%s.%s.log', source, mfilename);
if exist([source '/' fileOutLog],'file'),    delete([source '/' fileOutLog]);    end
if exist(source,'dir'),                       diary([source '/' fileOutLog]);    
else                        error(sprintf('Directory %s not found!',source));    end
% CHECK that the xml file is present
if ~exist([xml],'file'),    error('Parameter file %s not found!', xml);          end

% LOAD processing parameters from a xml file
par = LoadXml([xml]);
subSessionList = par.SubsessionDescription;

%NOTE: spike processing will be modified to switch to Kilosort
%Process spikes only if spike groups are defined in the xml file

prcSpikeFlag = isfield(par, 'SpkGrps') &&  ~isempty(par.SpkGrps(1).Channels);


% RENAME data files to conform to the ndm_naming scheme
% MOVE data files from subsession directories to the same session directory
%! ecube_link_files
%! ecube_rename_files

%%%>>>

%%%<<< Processing using MATLAB

% CONVERT open ecube .continous files to a .dat 
dat = fullfile(filebase,[filebase,'.dat']);
if ~exist(dat,'file') || overwrite,
    ecube_map_continuous_to_dat(filebase, xml, acqSystem, chanmap, processorList, subSessionList);
else
    fprintf('File %s found. Function oephys2dat was SKIPPED.\n', dat)
end

% COPY xml template to raw and processed directories
copyfile(xml, fullfile(filebase,[filebase,'.xml']), 'f');
%copyfile(xml, fullfile(source,  [filebase,'.xml']), 'f');
  
% COPY channel map template to raw and processed directories
if ~isempty(channelMap)
    copyfile(chanmap, fullfile(filebase, [filebase, '.chanmap']), 'f');
    copyfile(chanmap, fullfile(source,   [filebase, '.chanmap']), 'f');    
end

%%%>>>



%%%<<<  MOVE processed files to final 'processed' folder
% CLEAN up files
%CleanProcessedSession(filebase)    
cwd = pwd();
finalDirPath = fullfile('..','..','..','processed','ephys',subject,filebase);
system(['mv ',filebase,' ',finalDirPath]);
cd(fullfile(finalDirPath))

%%%>>>


%%%<<< Processing using NDM scripts
% CREATE a LFP file by downsampling the .dat file
lfp = fullfile(filebase,[filebase,'.lfp']);
if ~exist(lfp,'file')
    [status, cmdout] = system(sprintf('ndm_lfp %s.xml', filebase), '-echo');
else
    fprintf('File %s found. Function ndm_lfp was SKIPPED.\n', lfp)
end

% PROCESS spikes only if necessary (NOT TESTED YET)
if prcSpikeFlag==1
% CREATE a .fil file by highpass filtering the .dat file (for future spike detection)
    [status, cmdout] = system(sprintf('ndm_hipass %s.xml', filebase), '-echo');
% DETECT spikes in the .fil file (requires defined spike groups in the xml file!)
    [status, cmdout] = system(sprintf('ndm_extractspikes %s.xml', filebase), '-echo');
% COMPUTE PCA for the detected spikes
    [status, cmdout] = system(sprintf('ndm_pca %s.xml', filebase), '-echo');

    for s = 1:numel(par.SpkGrps),
        fets = strrep(num2str(ones([1,par.SpkGrps(s).nFeatures...
                                      *numel(par.SpkGrps(s).Channels)+1])),' ','');
        system(['nohup KlustaKwik ',filebase,' ',num2str(s),                  ...
                ' -MinClusters 50 -MaxClusters 200 -MaxPossibleClusters 200 ' ...
                '-UseFeatures ' fets ' &']);
    end
end %if prcSpikeFlag==1

cd(cwd); 

%%%>>>

diary off
% MOVE ProcessOephys log file to the directory with converted files
if exist(fullfile(source,fileOutLog), 'file')
    movefile(fullfile(source,fileOutLog), fullfile(finalDirPath,fileOutLog));
end
