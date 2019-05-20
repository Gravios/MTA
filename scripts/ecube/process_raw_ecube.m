function process_raw_ecube(Session, xml, varargin)
% function  process_raw_ecube(filebase, xml, <chanMapFile>, <processorList>)
%
% PREPARE data for processing:
%     Place all subsessions in subject directory { i.e. /*/*/data/raw/ecube/IF13 }
%    
%    PLACE Place all subsessions in subject directory in /storage/<User>/data/raw/<AnimalID>/.
%    CREATE an ASCII *.note file with experiment's metadata. {template:  oephys/template_files/ }
%    PLACE *.note in processing directory. 
%    RENAME *.note template copy as AnimalID-YYYYMMDD.note 
%    UPDATE *.note with experiment's metadata
% 
%    CREATE XML parameter file for data processing. (AKA "parameter file")   { template: oephys/template_files/  }
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
% Dependencies:
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
% USAGE : process_raw_ecube(filebase, xml, <chanMapFile>, <processorList>)
%
%    INPUT:
%     Session     is a name of the open ephys session to be processed (AnimalID-YYYMMDD-oephys).
%     xml             is a XML file with parameters necessary for processing (such as subsession descriptions).
% <chanMapFile>       is an ASCII file with a list of open ephys channel indicies in a proper anatomical order (raw- or columnwise).
%                     if not present, channel mapping is skipped.
% <processorList>   is an ID of the open ephys recording processor from which the data files must be merged. Default=[];
%                     (e.g. a file 100_CH1.continuous was recorded by the processor 100).
%                     This input parameter can be skipped if data files from only one processor are present.
%                     If data files from different processors are present and processorList is not provided, ProcessOephys 
%                     will stop and show an error mesage.
%OUTPUT:
% The script creates a new directory AnimalID-YYYMMDD next to the directory with the raw data fiels AnimalID-YYYMMDD-oephys
% and saves there all the converted files: AnimalID-YYYMMDD.dat, AnimalID-YYYMMDD.lfp, AnimalID-YYYMMDD.cat.evt and additional 
% ASCII .info files with info about processed subsessions).
% All processing steps are logged into an ASCII file AnimalID-YYYMMDD-oephys.ProcessOephys.log.
%
%EXAMPLE:   ProcessEcube('TestEcube-20180619-ecube', '1ses-160chan.xml', [])
%
%DEPENDENCIES:
% labbox: DefaultArgs, LoadXml 
% oephys: CorrectNames, oephys2dat, CleanProcessedSession
% NDM toolbox (not matlab): ndm_lfp, ndm_hipass, ndm_extractspikes, ndm_pca
% 
% Evgeny Resnik
% version 10.09.2017
%
% 
%TO DO:
% - fix incorrect renaming of non-default files in the session folder (e.g. avi).
% - oephys2dat_subses_blocks - replace cell2mat(raw_t) with a loop as in oephys2dat_subses. 
% - oephys2dat_subses_blocks - make proper calculation of start/end times for time gaps if present in the data. 
% - convert open ephys event file
% - processing of bonsai video files





% Session = 'TestEcube-20180618-ecube'PowerPcnt
% xml = '4ses-96chan.xml';

% Session = 'TestEcube-20180619-ecube'
% xml = '1ses-160chan.xml';

% Session = 'TestEcube-20180620-ecube'
% xml = '3ses-160chan.xml';

addpath('/storage/gerrit/code/matlab/Wraps/Convert/openephys/')
addpath('/storage/gerrit/code/matlab/Wraps/Convert/Ecube/')

Session     = 'IF13-20190410a-ecube';
xml     = 'IF13.xml';
chanMapFile = 'IF13.chanmap';
processorList = [100,103]; 
AcqSystem = 'ecube';

% chanMapFile = [];
% processorList = [];
% mfilename = 'ProcessEcube';



if nargin<1
    error(['USAGE:  ProcessEcube(Session, xml, <chanMapFile>, <processorList>)'])
end


%Parse input parameters
[ chanMapFile, processorList] = DefaultArgs(varargin,{ [], [] });


fprintf(' \n')
fprintf(' \n')
fprintf('=======================================================================================================\n')
fprintf(['                             %s - %s\n'], Session, mfilename)
fprintf('=======================================================================================================\n')



AcqSystem = 'ecube';


%----------------------------------------------------------------------------------------------------%
%                                    Preparations
%----------------------------------------------------------------------------------------------------%

%Extract filebase from the session+prefix name
out = regexp(Session, '-', 'split');
% % if length(out)<3 || ~strcmp(out{3}, 'ecube')
% %     error('Session directory must be named as "AnimalID-YYYMMDD-ecube".')
% % end
fileBase = [out{1} '-' out{2}];
clear out

%Check that the function is run outside the session directory
out = regexp(pwd, '/','split'); 
if any(ismember(out, Session)) || any(ismember(out, filebase)) 
    cd('..')
    fprintf('This function must be run outside the directory to be processed!\n')
end
clear out

%make new directory named by filebase if it doesn't exist already
%mkdir(filebase);

%Start logging the command window messages
fileOutLog = sprintf('%s.%s.log', Session, mfilename);
if exist([Session '/' fileOutLog],'file')
   delete([Session '/' fileOutLog])
end
if exist(Session,'dir')
    diary([Session '/' fileOutLog])
else
    error(sprintf('Directory %s not found!',Session))
end


%Check that the xml file is present
if ~exist([xml],'file')
    error('Parameter file %s not found!', xml)
end


%Load processing parameters from a xml file
%par = LoadXml([filebase '/' xml]);
par = LoadXml([xml]);
SubSesDescriptions = par.SubsessionDescription;



%NOTE: spike processing will be modified to switch to Kilosort
%Process spikes only if spike groups are defined in the xml file
if isfield(par, 'SpkGrps') &&  ~isempty(SpikeGroups(1).Channels)
    SpikeGroups = par.SpkGrps;
    prcSpikeFlag=1;
else
    SpikeGroups = [];
    prcSpikeFlag=0;
end


%----------------------------------------------------------------------------------------------------%
%                                   Processing using MATLAB
%----------------------------------------------------------------------------------------------------%

% RENAME data files to conform to the ndm_naming scheme
% MOVE data files from subsession directories to the same session directory
%! ecube_link_files
%! ecube_rename_files
% CorrectNames(Session, AcqSystem, SubSesDescriptions)
 

%Convert open ephys .continous files to a .dat file for individual subsessions
%Merge .dat files from the subsessions to a single file
file = sprintf('%s/%s.dat', filebase, filebase);
if ~exist(file,'file')
    map_oephys_to_dat(Session, xml, AcqSystem, chanMapFile, [], processorList)
else
    fprintf('File %s found. Function oephys2dat was SKIPPED.\n', file)
end




%Copy the template xml file to the directories with raw and processed data files
file = sprintf('%s/%s.xml', filebase, filebase);
if ~exist(file,'file')
    copyfile(xml, file, 'f')
else
    fprintf('File %s found. Copying %s --> %s was SKIPPED.\n', file, xml, file)
end

file = sprintf('%s/%s.xml', Session, filebase);
if ~exist(file,'file')
    copyfile(xml, file, 'f')
else
    fprintf('File %s found. Copying %s --> %s was SKIPPED.\n', file, xml, file)
end
  


%Copy the template channel mapping file to the directories with raw and processed data files
if ~isempty(chanMapFile)
    
    [~, ~, ext] = fileparts(chanMapFile);   
   
    file = sprintf('%s/%s%s', filebase, filebase, ext);
    if ~exist(file,'file')
        copyfile(chanMapFile, file, 'f')
    else
        fprintf('File %s found. Copying %s --> %s was SKIPPED.\n', file, chanMapFile, file)
    end
    
    file = sprintf('%s/%s%s', Session, filebase, ext);
    if ~exist(file,'file')
        copyfile(chanMapFile, file, 'f')
    else
        fprintf('File %s found. Copying %s --> %s was SKIPPED.\n', file, chanMapFile, file)
    end   
    
end



% % return



%----------------------------------------------------------------------------------------------------%
%                                   Processing using NDM scripts
%----------------------------------------------------------------------------------------------------%
%Enter the directory with processed data files
cd(filebase)

%Create a LFP file by downsampling the .dat file
file = sprintf('%s.lfp', filebase);
if ~exist(file,'file')
    [status, cmdout] = system(sprintf('ndm_lfp %s.xml', filebase), '-echo');
else
    fprintf('File %s found. Function ndm_lfp was SKIPPED.\n', file)
end


% PROCESS spikes only if necessary (NOT TESTED YET)
if prcSpikeFlag==1
% CREATE a .fil file by highpass filtering the .dat file (for future spike detection)
    [status, cmdout] = system(sprintf('ndm_hipass %s.xml', filebase), '-echo');
% DETECT spikes in the .fil file (requires defined spike groups in the xml file!)
    [status, cmdout] = system(sprintf('ndm_extractspikes %s.xml', filebase), '-echo');
% COMPUTE PCA for the detected spikes
    [status, cmdout] = system(sprintf('ndm_pca %s.xml', filebase), '-echo');
end %if prcSpikeFlag==1


%Exit the directory with processed data files
cd('..')



%----------------------------------------------------------------------------------------------------%
%            Deleting unnecessary files from the directory with processed data files
%----------------------------------------------------------------------------------------------------%
CleanProcessedSession(filebase)    

%move processed files to final 'processed' folder

out = regexp(pwd, '/', 'split');
new = [];
for k=1:length(out);
    if ~strcmp(out{k},'');
        if ~strcmp(out{k},'raw');
        new=[new '/' out{k}];
        else
            new=[new '/processed'];
        end       
    end
end

newdir=[new '/' filebase];
clear out new
mkdir2(newdir);

system(sprintf('mv %s/* %s', filebase, newdir), '-echo')

system(sprintf('rm -r  %s/', filebase), '-echo')
system(sprintf('rm %s.chanmap', filebase), '-echo')
system(sprintf('rm %s.xml', filebase), '-echo')

%Exit the directory with processed data files
cd('..')



%---------------------------------------- END ---------------------------------------------------------%
diary off

%Move ProcessOephys log file to the directory with converted files
if exist([Session '/' fileOutLog], 'file')
    movefile([Session '/' fileOutLog], [filebase '/' fileOutLog])
end



return

% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %------------------------------------------------------------------------------------------------------------%
% $$$ %                     Interactive steps after running ProcessOephys on raw open ephys data
% $$$ %                    (for details see the lab wiki "Initial processing of open ephys data")
% $$$ %------------------------------------------------------------------------------------------------------------%
% $$$ % -Move the directory with processed data files to /storage/<User>/data/processed/<AnimalID>/
% $$$ 
% $$$ % -Check filebase.dat and filebase.lfp files in Neuroscope:
% $$$ %    - mark permanently bad channels as 'skipped'
% $$$ %    - identify CA1pyr channel if present
% $$$ 
% $$$ % 
% $$$ % -Create an ASCII .mypar file with all the parameters necessary for further data processing.
% $$$ %  For this copy a template .mypar file from oephys/template_files/ to the directory with processed data files. 
% $$$ %  Rename the file as AnimalID-YYYYMMDD.mypar and update it with actual information using any text editor (e.g. gedit, emacs). 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %----------------------------------------------------------------------------------------------------%
% $$$ %                        Common processing steps for any experimental session
% $$$ %----------------------------------------------------------------------------------------------------%
% $$$ %Enter the directory with processed data files
% $$$ cd(filebase)
% $$$ 
% $$$ %Create a .lfpinterp file with bad channels (marked as skipped in .xml file) replaced with interpolated values
% $$$ %Calculate the mean coherence between channels to visualize possible cross-talks
% $$$ %Calculate time maps of spectral similarity between channels to detect possible bad periods
% $$$ ProcessClean(filebase)
% $$$ 
% $$$ %Create .sts files with start/end timestamps of run/sleep behavioral episodes pooled together by state
% $$$ %I should rename it to ProcessBhvEpisodes or inlcude to ProcessClean
% $$$ % ProcessStates(CurrentFileBase)
% $$$ 
% $$$ 
% $$$ % mypar is required starting from here!!
% $$$ 
% $$$ %Detect brain states (RUN, SWS, REM) using theta power and, if available, speed or accelerometer data
% $$$ DetectBrainStates(filebase)
% $$$ 
% $$$ %Visual check-up of the detected brain state periods
% $$$ CheckEegStates(filebase)
% $$$ 
% $$$ %Plot a figure with the final brain state periods
% $$$ PlotBrainStates(filebase)
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ %------------------------------------------------------------------------------------------------------------%
% $$$ %                                   Detect ripples, compute parameters of ripples and theta                
% $$$ %------------------------------------------------------------------------------------------------------------%
% $$$ %I am going to split theta and ripples processing parts into different functions
% $$$ %Detect theta and ripples
% $$$ ProcessThetaRipples(CurrentFileBase)
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
