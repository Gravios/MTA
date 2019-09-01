function ecube_map_continuous_to_dat(filebase, xml, varargin)
% function map_oephys_to_dat(Session, xml, acqSystem, channelMap, processorList, subSessionList)
%
% CONVERTS .continuous files with wideband signal recorded 
%     by an open ephys acquisition acqSystem to .dat data format.
%     NOTE : This function must be applied to an experimental session directory
%     which contains data files that have already been renamed to conform the ndm
%     naming scheme by CorrectNames.m.
%
% INPUT :
%     filebase - string: name of the session to be processed {format AnimalID-YYYMMDD-oephys}.
%     xml -      string: xml file name, file contains parameters necessary for processing
%     acqSystem -   string: name of the aqcuisition system 
%                        {'oephys', 'openephys'; 'ecube', 'eCube'; 'nlx', 'Neuralynx'}
%     channelMap, string: Name of ASCII file, which contains list of channel indicies in 
%                 anatomical order (raw- or columnwise).
%                 IF not present, channel mapping is skipped.
%     processorList -  cellArray: list of open ephys pipeline processor IDs to be merged.
%     subSessionList - cellArray: list of subsessions whose .dat files must be merged. 
%
% USAGE : 
%     oephys2dat(filebase, xml, acqSystem, channelMap, processorList, subSessionList)
%
% OUTPUT :
%     savefile - AnimalID-YYYYMMDD.dat 
%
% EXAMPLE :    
%     oephys2dat('ER50-20170101-oephys', 'template-A32.xml', 'oephys')
%
% DEPENDENCIES :
%     labbox: DefaultArgs, memorylinux, mkdir2
%     oephys: load_open_ephys_header, oephys2dat_subses
%     FMAToolbox: SaveBinary
%
% Evgeny Resnik
% version 04.12.2017
% version 07.06.2018: can deal now with eCube data
% Justin Graboski
% version 29.07.2019: refactoring
%
%--------------------------------------------------------------------------------------------------%
%                               OBSERVED ISSUES WITH eCube/OpenEphys DATA:
% 1) Normally successive timestamps within the same .continuous file increment by one sample.
%    In some case however, there might be gap(s). Importantly, these gaps can differ between 
%    CH and PAI channels. If the data contain time gaps (>1 .dat sample), a new time vector
%    without gaps is created.
% 
% 2) Data recorded by different processors (pipeline modules in the software, like CH, PAI, PDI
%    in eCube or CH, AUX - in OpenEphys) may have different number of samples. In other cases 
%    they have same sample length, but shifted (by 10-60 ms) timestamps.
% 
%--------------------------------------------------------------------------------------------------%


% DEFARGS ------------------------------------------------------------------------------------------

assert(nargin>2,'USAGE:  oephys2dat(Session, xml, channelMap, subSessionList, acqSystem)')

[ channelMap, subSessionList, acqSystem ] = DefaultArgs(varargin,{  [], [], 'ecube' });
%---------------------------------------------------------------------------------------------------


%%%<<< TESTING VARS
% $$$ addpath('/storage/gerrit/code/matlab/Wraps/Convert/openephys/')
% $$$ addpath('/storage/gerrit/code/matlab/Wraps/Convert/Ecube/')
% $$$ 
% $$$ filebase       = 'IF13-20190409c';
% $$$ xml            = 'IF13.xml';
% $$$ channelMap     = 'IF13.chanmap';
% $$$ processorList  = {'100','103'}; 
% $$$ acqSystem      = 'ecube';
% $$$ mfilename      = 'map_oephys_to_dat';
% $$$ subSessionList = 'ALL';
%%%>>>


% SET constant parameters
% Maximum size of a sesion the PC can deal with (used in the check on size of future .dat files)
% Based on my empirical testing, our server can not load all .conitnuous files from a given
% subsession to RAM, if their overall size exceeds this value. In this case the files are 
% converted by smaller time chunks (see below).
% MaxDatFileSize_Gb = 20;
MAX_DATFILE_SIZE = 15;
source = [filebase,'-',acqSystem];
% ASSERT xml file existance
assert(exist(xml,'file')==2,sprintf('Parameter file %s not found', xml))


fprintf(['%s:  %s (converting and merging .continous --> .dat) \n'], filebase, mfilename );
fprintf(['%s \n'], datestr(clock) );


%%%<<< GET processor list

% NOT necessary since processors are included in channelMap
% $$$ sys = dir(fullfile(source,[filebase,'*-settings.xml']));
% $$$ sys = xml2struct(fullfile(source,sys(1).name));
% $$$ assert(~isempty(sys),'map_oephys_to_dat:XMLNotFound');
% $$$ processorList = ecube_select_ephys_processors(sys);

%%%>>>

%%%<<< LOAD channel mapping from the ASCII file    

if exist(channelMap, 'file'),
    fprintf('Loading channel mapping from %s ->', channelMap);
    fid = fopen(channelMap);
    chanmap = textscan(fid,'%s');
    fclose(fid);    
    chanmap = chanmap{1};
    fprintf(' DONE\n');
    if isempty(chanmap),
        fprintf('Channel Map vector is empty -> MAPPING all channels in original order..\n')    
    end
else
% Could also parse the channel order from the settis.xml 
% settings.signalchain.processor: 100 channels
%                                 101 channel map
%                                 102 lfp viewer
%                                 103 analogue input
%    <PROCESSOR name="Filters/Channel Map" insertionPoint="1" pluginName="Channel Map"
%               pluginType="1" pluginIndex="3" libraryName="Channel Mapper" libraryVersion="1"
%               isSource="0" isSink="0" NodeId="101">
%        <EDITOR isCollapsed="0" displayName="Channel Map" Type="ChannelMappingEditor">
    fprintf('Channel Map file is absent -> MAPPING all channels in original order. \n')
    chanmap = [];
end

%%%>>>

%%%<<< CREATE a list of all .continious files

dirList = dir(fullfile(source,[filebase,'*.continuous']));
filenames = {dirList.name}';
assert(~isempty(filenames),'map_oephys_to_dat:EphysFilesNotFound');
clear('dirList');

%%%>>>

%%%<<< PARSE file names to extract information about the files

selectedChannels = false(size(filenames));
chanmapInd = [];
clear('chanInfo');
for k=1:length(filenames)
    cit= regexp(filenames{k},['(?<subjectName>[a-zA-Z]{1,4}\d+)-',...
                              '(?<sessionName>\d{8,8}[a-zA-Z]?)-',...
                              '(?<subSesId>\d+)-',                ...
                              '(?<subSesName>[a-zA-Z]+)-',        ...
                              '(?<acqSystem>\w+)-',               ...
                              '(?<processorId>\d\d\d+)-',         ...
                              '(?<chanType>\w+)-',                ...
                              '(?<chanId>\w+)\.',                ...
                              '(?<ext>.*)'],'names');
    
    cit.chanName = [cit.processorId,'-',cit.chanType,'-',cit.chanId];
    
    if ~isempty(chanmap)
        cmi = find(~cellfun(@isempty,regexp(cit.chanName,cf(@(c) ['^',c,'$'],chanmap),'match')));
    else
        cmi = k;
    end

    if isempty(cmi)
        continue
    end
    
    chanInfo(cmi) = cit;
    chanmapInd(end+1) = cmi;        
    selectedChannels(k) = true;
    
end
% $$$ clear('out','out2','k');
clear('cit','k');

assert(numel(unique(chanmapInd))==numel(chanmap),...
       'map_oephys_to_dat:chanmap:DuplicateChannel');


subSesIdList = unique({chanInfo.subSesId});
nSubSes = length(subSesIdList);
subSesNameList = unique({chanInfo.subSesName});
chanTypeList = unique({chanInfo.chanType});

selectedFiles = reshape(filenames(selectedChannels),numel(chanmap),nSubSes);
for s = 1:nSubSes,
% REORDER files in accordance to chanmap
    selectedFiles(chanmapInd,s) = selectedFiles(:,s);
end

%%%>>>

%%%<<< CHECK on size of future .dat files

% CHECK if size of individual subsesions data exceeds server limit

% COMPUTE Size (Gb) of individual sibsessions
for s=1:nSubSes
% SIZE (Gb) of all files from the given subsession
    if isunix && ~unix(['test -L ',fullfile(source,selectedFiles{1,s})])
        [~,out] = unix(['readlink ',fullfile(source,selectedFiles{1,s})]);
        out(double(out)<=32) = [];        
        out = dir(out);
    else
        out = dir(fullfile(source,selectedFiles{1,s}));    
    end
    subSesSizeGb(s) = out.bytes/1e9 * size(selectedFiles,1);
    clear out 
end 

%Size (Gb) of all the subsessions
subSesSizeGbAll = sum(subSesSizeGb);

% CHECK whether any of the subsession exceeds the maximum size
if any(subSesSizeGb > MAX_DATFILE_SIZE)  || subSesSizeGbAll > MAX_DATFILE_SIZE  
    for s=1:nSubSes
        fprintf('%s: %1.1f Gb \n', subSesNameList{s}, subSesSizeGb(s));
    end
    fprintf('COMPUTING future .dat file size ... %1.1f Gb \n', subSesSizeGbAll);
end

% ??? For now it doesn't deal with multiple files with suffixes because of timestamp issues
% $$$ if ~all( cellfun(@(x) isempty(x), postfix) )
% $$$     error('ERROR: found multiple files with suffixes for each channel. Stopped. Contact Evgeny')
% $$$ end

%%%>>>

%%%<<< PROCESS individual subsessions 
% CREATE an output directory for converted files
create_directory(filebase);

cwd = pwd();
cd(source)

%Loop across subsessions
for s=1:nSubSes    
    
    % name of the output .dat file
    if nSubSes==1,
        dat = sprintf(['%s.dat'], filebase);
    else
        dat = sprintf(['%s-%02.0f-%s.dat'], filebase, subSesIdList(s), subSesNameList{s});
    end
    
% MAINTAIN list of individual .dat files for future merging
    datList{s} = dat;    

    fprintf('CONVERTING .continuous files -> %s \n', dat)
    
% SKIP if the file exists
    if exist(dat, 'file')
        fprintf('File %s already exists. Converting was skipped. \n', dat)
        
        %move the .dat and .dat.sts files to the dedicated directory
        if exist(dat, 'file')
            movefile(dat, fullfile(cwd, filebase, dat));
        end
        if exist([dat '.info'], 'file')
            movefile([dat '.info'], fullfile(cwd, filebase, [dat '.info']));
        end 
        if exist([dat '.ChunksByGaps'], 'file')
            movefile([dat '.ChunksByGaps'], fullfile(cwd, filebase, [dat '.ChunksByGaps']));
        end
        
        continue
        
    elseif exist(fullfile(cwd,filebase,dat), 'file'),
        fprintf('File %s already exists. Converting was skipped. \n', ...
                fullfile(cwd, filebase, dat));
        continue
    end
        
    ecube_map_continuous_to_dat_subses_block( filebase, selectedFiles(:,s), dat, chanInfo, acqSystem, MAX_DATFILE_SIZE);
    
    %move the .dat and .dat.sts files to the dedicated directory
    if exist(dat, 'file')
        movefile(dat,fullfile(cwd,filebase, dat))
    end    
    if exist([dat '.info'], 'file')
        movefile([dat '.info'], fullfile(cwd, filebase, [dat '.info']))
    end
    if exist([dat '.ChunksByGaps'], 'file')
        movefile([dat '.ChunksByGaps'], fullfile(cwd, filebase, [dat '.ChunksByGaps']))
    end
    
    clear('dat');
end%for s

cd(cwd);

%%%>>>

%%%<<< MERGE .dat files across subsessions 
%Merge .dat files from individual subsessions
%Create filebase.cat.evt with start/stop timestamps (ms) of individual subsessions
if nSubSes > 1,
    cd(filebase)
    fileMergedOut = sprintf('%s.dat', filebase);
% SKIP if the file exists
    if ~exist(fileMergedOut, 'file')
        fprintf(' mergedat: merging subsession .dat files ---> %s \n', fileMergedOut)        
        tic
        mergedat(datList, fileMergedOut)
        toc
        fprintf(' DONE\n');        
    else
        fprintf('File %s already exists. Merging was skipped. \n', fileMergedOut)
    end
    cd('..')    
end

%%%>>>
diary off





    
    



