function map_oephys_to_dat(filebase, XmlFile, AcqSystem, varargin)
%oephys2dat is a function which converts .continuous files with wideband signal recorded 
%by an open ephys acquisition system to .dat data format.
%NOTE: This function must be applied to an experimental session directory
%which contains data files that have already been renamed to conform the ndm
%naming scheme by CorrectNames.m.
%
%USAGE: oephys2dat(filebase, XmlFile, AcqSystem, <ChanMapFile>, <SubSes2Merge>, <Processor2Merge>)
%
%INPUT:
% filebase     is a name of the session to be processed in the following format AnimalID-YYYMMDD-oephys.
% XmlFile             is a xml file with parameters necessary for processing such as subsession description 
%                     and channel mapping.
% AcqSystem           is a name of the aqcuisition system ('oephys' - open ephys; 'ecube' - eCube; 'nlx' - Neuralynx).
% <SubSes2Merge>      is a cell vector with a list of subsessions whose .dat files must be merged. Default=[] (all subsessions).
%                     Use it when some of the subsessions must be excluded from merging, i.e. because of too large files. 
%                     For example, a long (hours) sleep session should be processed individually as a separate session.
% <Processor2Merge>   is a open ephys pipeline processor ID whos data files must be merged. Default-the one found in the data files.
%
%OUTPUT:
% AnimalID-YYYMMDD.dat (mkV)
%
%
%EXAMPLE:    oephys2dat('ER50-20170101-oephys', 'template-A32.xml', 'oephys')
%
%DEPENDENCIES:
% labbox: DefaultArgs, memorylinux, start_parallel, mkdir2
% oephys: load_open_ephys_header, oephys2dat_subses
% FMAToolbox: SaveBinary
%
% Evgeny Resnik
% version 04.12.2017
% version 07.06.2018: can deal now with eCube data
%
%-------------------------------------------------------------------------------------------------------------------------------%
%                               OBSERVED ISSUES WITH eCube/OpenEphys DATA:
% 1) Normally successive timestamps within the same .continuous file increment by one sample. In some case however, 
%    there might be gap(s). Importantly, these gaps can differ between CH and PAI channels. If the data contain time gaps
%    (>1 .dat sample), a new time vector without gaps is created.
% 
% 2) Data recorded by different processors (pipeline modules in the software, like CH, PAI, PDI in eCube or CH, AUX - in OpenEphys) 
%    may have different number of samples. In other cases they have same sample length, but shifted (by 10-60 ms) timestamps.
% 
%-------------------------------------------------------------------------------------------------------------------------------%
% 
%TO DO:
% - CorrectNames: check whether subsession must be ordered based on their creation time or time extracted from the directory names.
% - If data files are in the main session directory treat it as a single  subsession experiment





% filebase = 'Exam4-20180606-ecube'
% Processor2Merge = [];
% SubSes2Merge = [];
% AcqSystem = 'ecube';
% ChannelMapping = [];
% XmlFile = 'ses3.xml';
% mfilename = 'oephys2dat';





if nargin<1
    error(['USAGE:  oephys2dat(filebase, XmlFile, AcqSystem, <ChanMapFile>, <SubSes2Merge>, <Processor2Merge>)'])
end


% Parse input parameters
[ChanMapFile, SubSes2Merge, Processor2Merge] = DefaultArgs(varargin,{ [], [], [] });


%---------------------------- Constant parameters ------------------------------------------------%
%Maximum size of a sesion the PC can deal with (used in the check on size of future .dat files)
%Based on my empirical testing, our lab server can not load all .conitnuous files from a given subsession to RAM,
%if their overall size exceeds this value. In this case the files are converted by smaller time chunks (see below).
% MaxDatFileSize_Gb = 20;
MaxDatFileSize_Gb = 15;

%-------------------------------------------------------------------------------------------------%
%Extract FileBase from the session name
out = regexp(filebase, '-', 'split');
if length(out)<3
    error('Session directory must be named as "AnimalID-YYYMMDD-<AcqSystem>".')
end
if ~ismember(out{3}, {'ecube', 'oephys','nlx'} )
    error('Session directory must be named as "AnimalID-YYYMMDD-<AcqSystem>, where AcqSystem must be either "oephys", "ecube" or "nlx".')
end
FileBase = [out{1} '-' out{2}];
clear out


%------------------------------ Checks ------------------------------------------------------%
%Check that the function runs outside the session directory
out = regexp(pwd, '/','split'); 
if any(ismember(out, filebase)) || any(ismember(out, FileBase))    
    error('This function must be run outside the directory to be processed!')
end
clear out


%Check that the xml file is present
if ~exist(XmlFile,'file')
    error(sprintf('Parameter file %s not found!', XmlFile))
end
    

if ~ismember(AcqSystem, {'oephys', 'ecube', 'nlx'})
    error(['Unknown acquisition system "' AcqSystem '"! Must be either "oephys", "ecube" or "nlx".'])
end



% %start logging the command window messages
% FileOutLog = sprintf('%s.%s.log', filebase, mfilename);
% if exist([filebase '/' FileOutLog],'file')
%    delete([filebase '/' FileOutLog])
% end
% diary([filebase '/' FileOutLog])




fprintf(' \n')
fprintf(' \n')
fprintf('=======================================================================================================\n')
fprintf('            %s:  %s (converting and merging .continous --> .dat) \n', filebase, mfilename )
fprintf('=======================================================================================================\n')
fprintf('                                 %s \n', datestr(clock) )
fprintf('=======================================================================================================\n')



% %OPTION-1: Load channel mapping from a XML file (not handy when the channel count is high)
% fprintf('Loading channel mapping from %s ...', XmlFile)
% par = LoadXml_ER(XmlFile);
% ChannelMapping = par.ChannelMapping;
% fprintf('DONE\n')


%OPTION-2: Load channel mapping from the ASCII file
if exist(ChanMapFile, 'file')
    fprintf('Loading channel mapping from %s ...', ChanMapFile)
    ChannelMapping = dlmread(ChanMapFile);
    fprintf('DONE\n')
else
    fprintf('Channel mapping file is not provided. No channel remapping is done. \n')
    ChannelMapping = [];
end


%Ensure that the channel mapping vector is a row (not column)
ChannelMapping = ChannelMapping(:);
ChannelMapping = ChannelMapping';


if isempty(ChannelMapping)
    fprintf('Channel mapping vector is empty. No channel remaping is done!\n')
end

%This check works only when channel indices are sequential, which is not
%the case when only a subset of the available channels are kept.
% if length(ChannelMapping)~=max(ChannelMapping)
%     error(sprintf('Channel mapping vector seems to be wrong!\n'))
% end


%Create a list of all .continious files in the session directory
out = dir([filebase '/*.continuous']);
out([out(:).isdir]) = [];
FileNames = arrayfun(@(x)x.name, out,'UniformOutput',0);
if length(FileNames)<1
    disp('No .continuous files found.');
    return;
end
clear out


%Check whether the files have been renamed by CorrectNames.m
out = regexp(FileNames{1}, '-', 'split');
if length(out)==1
    error('Data files must be first renamed by the function CorrectNames!')
end
clear out



%Parse file names to extract information about the files
% out = cellfun(@regexp, FileNames, repmat({'-'}, size(FileNames)), repmat({'split'}, size(FileNames)), 'UniformOutput', false);
clear SubSesID SubSesDescription ProcessorID ChanType Chan Postfix
for k=1:length(FileNames)
    out = regexp(FileNames{k}, '-', 'split');
    SubSesID(k)          = str2num(out{3});
    SubSesDescription{k} = out{4};
    ProcessorID(k)       = str2num(out{6});
    ChanType{k}          = out{7};  
    %parse depending on presence of a postfix
    if length(out)==8
        %file name without a postfix
        out2 = regexp(out{8}, '[.]', 'split');
        Chan(k) = str2num(out2{1});
        Postfix{k} = [];
    elseif length(out)==9
        %file name with a postfix
        Chan(k) = str2num(out{8});
        out2 = regexp(out{9}, '[.]', 'split');
        Postfix{k} = out2{1};
    end
end
clear out out2 k


%Lists of unique values
List_SubSesID = unique(SubSesID);
nSubSes = length(List_SubSesID);
List_ProcessorID = unique(ProcessorID);
List_ChanType = unique(ChanType);
for k=1:nSubSes
    ind = find(SubSesID == List_SubSesID(k), 1);
    List_SubSesDescription{k} = SubSesDescription{ind};
    clear ind
end




if strcmp(AcqSystem, 'oephys')
    %OpenEphys system in our lab usually saves data files from only one pipeline module (processor),
    %therefore, generate an error message if multiple processors are present.
    
    %Determine which processor must be used
    if isempty(Processor2Merge) && length(List_ProcessorID)==1
        %if processor is not provided as an input parameter and only one processor is found in the data, use that one
        Processor2Merge = List_ProcessorID;
    elseif  isempty(Processor2Merge) && length(List_ProcessorID)>1
        %if processor is not provided as an input parameter and multiple processors are found in the data
        error('Data files recorded by multile OpenEphys processors found. Please specify which processor must be used in Processor2Merge!')
    end    
    
    %Discard files recorded by the "wrong" processors from the processing list
    if  length(List_ProcessorID)>1 && all(ismember(Processor2Merge, List_ProcessorID))
        fprintf('WARNING: .continuous files recorded by multiple open ephys processors found! Processing only files from processor "%d".\n', Processor2Merge)
        BadFiles = ProcessorID~=Processor2Merge;
        ProcessorID(BadFiles)       = [];
        SubSesID(BadFiles)          = [];
        SubSesDescription(BadFiles) = [];
        ChanType(BadFiles)          = [];
        Chan(BadFiles)              = [];
        Postfix(BadFiles)           = [];
        FileNames(BadFiles)         = [];
        clear BadFiles
    elseif ~all(ismember(Processor2Merge, List_ProcessorID))
        error( sprintf('.continuous files recorded by open ephys processor "%d" not found!', Processor2Merge) )
    end
    
        
elseif strcmp(AcqSystem, 'ecube')
    %eCube system, in contrast to OpenEphys, saves headstage channels, analog input and digital input channels via different pipeline modules
    %(processors), therefore merge dats files from all detected processors together
    disp('')
    fprintf('eCube system saves headstage channels, analog and digital input channels via different pipeline processors.\n')
    fprintf('Therefore, the parameter Processor2Merge is ignored and data from all found processors are merged together.\n')
    disp('')
    
    %This parameter will be further ignored in case of eCube file
    Processor2Merge = [];     
    
    %TO DO: check that channels of the same type (CH, PAI, PDI) have then same processor!
end %if strcmp(AcqSystem, 'oephys')



%Discard files belonging to the subsessions to be excluded from the processing list
if ~isempty(SubSes2Merge)
    
    %check that there is no repetitions
    if length(unique(SubSes2Merge)) ~= length(SubSes2Merge)
        error('SubSes2Merge contains repeating same subsession descriptions!')
    end
    
    BadSubSes = find(ismember(SubSes2Merge, List_SubSesDescription)==0);
    if ~isempty(SubSes2Merge) && isempty(BadSubSes)
        SubSes2Discard = setxor(List_SubSesDescription, SubSes2Merge);
        fprintf('WARNING: .continuous files from subsession "%s" are NOT processed (SubSes2Merge is not empty). \n', SubSes2Discard{:} )
        BadFiles = ismember(SubSesDescription, SubSes2Discard);
        ProcessorID(BadFiles)       = [];
        SubSesID(BadFiles)          = [];
        SubSesDescription(BadFiles) = [];
        ChanType(BadFiles)          = [];
        Chan(BadFiles)              = [];
        Postfix(BadFiles)           = [];
        FileNames(BadFiles)         = [];
        
        %re-compute lists of unique subsession descriptions
        List_SubSesID = unique(SubSesID);
        nSubSes = length(List_SubSesID);
        clear List_SubSesDescription
        for k=1:nSubSes
            ind = find(SubSesID == List_SubSesID(k), 1);
            List_SubSesDescription{k} = SubSesDescription{ind};
            clear ind
        end
        
    elseif ~isempty(BadSubSes)
        error( sprintf('.continuous files for this subsession not found: %s\n', SubSes2Merge{BadSubSes}) )
    end
    clear BadSubSes BadFiles SubSes2Discard
    
end %if ~isempty(SubSes2Merge)



%----------------------- Check on size of future .dat files ------------------------------------------------%
%Check wether the overall size of data in individual subsesions exceeds 40Gb (limits of our lab server nodes).

%Compute Size (Gb) of individual sibsessions
for s=1:nSubSes
    %list of subsession-specific .continuous files
    FileNames_Subses = FileNames(SubSesID==List_SubSesID(s));
    %Size (Gb) of all files from the given subsession
    out = dir([filebase '/' FileNames_Subses{1}]);
    SubSesSize_Gb(s) = out.bytes/1e9 * length(FileNames_Subses);
    clear FileNames_Subses out 
end 

%Size (Gb) of all the subsession
AllSubSesSize_Gb = sum(SubSesSize_Gb);

%Check whether any of the subsession exceeds the maximum size
if any(SubSesSize_Gb > MaxDatFileSize_Gb)  || AllSubSesSize_Gb > MaxDatFileSize_Gb  
    fprintf('-------------------------------------------------------------------------------------------------------\n')
    fprintf('Size of future .dat files in individual subsessions: \n');
    for s=1:nSubSes
        fprintf('%s: %1.1f Gb \n', List_SubSesDescription{s}, SubSesSize_Gb(s));
    end
    fprintf('Size of future .dat file merged across the subsessions: %1.1f Gb \n', AllSubSesSize_Gb);    
    fprintf('-------------------------------------------------------------------------------------------------------\n')    
end



%For now it doesn't deal with multiple files with suffixes because of timestamp issues
if ~all( cellfun(@(x) isempty(x), Postfix) )
    error('Multiple files with suffixes for each channel found. Stopped. Contact Evgeny!')
end


%----------------------- Process individual subsessions ------------------------------------------------%
%Create an output directory for converted files
mkdir2(FileBase)


cd(filebase)

%Loop across subsessions
for s=1:nSubSes    
    
    %list of subsession-specific .continuous files
    FileNames_Subses = FileNames(SubSesID==List_SubSesID(s));
   
    %name of the output .dat file
    FileOut = sprintf(['%s-%02.0f-%s.dat'], FileBase, List_SubSesID(s), List_SubSesDescription{s});
    
    %Keep list of individual .dat files for future merging
    Files2Merge{s} = FileOut;    
            
    fprintf('-------------------------------------------------------------------------------------------------------\n')
    fprintf('            oephys2dat_subses: converting .continuous files ---> %s \n', FileOut)
    fprintf('-------------------------------------------------------------------------------------------------------\n')
    
    %Skip if the file exists
    if exist(FileOut, 'file')
        fprintf('File %s already exists. Converting was skipped. \n', FileOut)
        
        %move the .dat and .dat.sts files to the dedicated directory
        if exist(FileOut, 'file')
            movefile(FileOut, ['../' FileBase '/' FileOut])
        end
        if exist([FileOut '.info'], 'file')
            movefile([FileOut '.info'], ['../' FileBase '/' FileOut '.info'])
        end 
        if exist([FileOut '.ChunksByGaps'], 'file')
            movefile([FileOut '.ChunksByGaps'], ['../' FileBase '/' FileOut '.ChunksByGaps'])
        end
        
        continue
        
    elseif exist(['../' FileBase '/' FileOut], 'file')
        %%exist(['../' FileBase '/' FileOut], 'file')
        
        fprintf('File %s already exists. Converting was skipped. \n', ['../' FileBase '/' FileOut])
        continue
    end
    
    
    %Convert .continuous to .dat for the given subsession files
    %-if the total size of all .continuous files <= MaxDatFileSize_Gb: load whole individual channels and them saave to the .dat file
    %-if the total size of all .continuous files >MaxDatFileSize_Gb: load and save to a .dat file in short blocks (1024 samples each)     
% %     if SubSesSize_Gb(s) <= MaxDatFileSize_Gb
% %         oephys2dat_subses(FileNames_Subses, FileOut, AcqSystem, ChannelMapping)
% %     else
% %         fprintf('The overall size of the subsession "%d-%s" (%1.0f Gb) exceeds the limit of %1.0f Gb. Chunk approach is used.\n', ...
% %             List_SubSesID(s), List_SubSesDescription{s}, SubSesSize_Gb(s), MaxDatFileSize_Gb)        
        oephys2dat_subses_blocks(FileNames_Subses, FileOut, AcqSystem, ChannelMapping, MaxDatFileSize_Gb) 
% %     end

    
% %TO BE DONE: use try-catch to catch memory errors here.
% if SubSesSize_Gb(s) <= MaxDatFileSize_Gb    
%     try
%         oephys2dat_subses(FileNames_Subses, FileOut, ChannelMapping)
%     catch ME
%         fprintf('ERROR: %s',ME.identifier)
%     end  
% else
%     fprintf('The overall size of the subsession "%d-%s" (%1.0f Gb) exceeds the limit of %1.0f Gb. Chunk approach is used.\n', ...
%         List_SubSesID(s), List_SubSesDescription{s}, SubSesSize_Gb(s), MaxDatFileSize_Gb)
%     oephys2dat_subses_blocks(FileNames_Subses, FileOut, ChannelMapping, MaxDatFileSize_Gb)
% end

    
    
    
    
    
    
    %move the .dat and .dat.sts files to the dedicated directory
    if exist(FileOut, 'file')
        movefile(FileOut, ['../' FileBase '/' FileOut])
    end    
    if exist([FileOut '.info'], 'file')
        movefile([FileOut '.info'], ['../' FileBase '/' FileOut '.info'])
    end
    if exist([FileOut '.ChunksByGaps'], 'file')
        movefile([FileOut '.ChunksByGaps'], ['../' FileBase '/' FileOut '.ChunksByGaps'])
    end
    
    clear FileNames_Subses FileOut 
end %loop across subsessions

cd('..')


% %Move oephys2dat log file to the directory with converted files
% if exist([filebase '/' FileOutLog], 'file')
%     movefile([filebase '/' FileOutLog], [FileBase '/' FileOutLog])
% end



%----------------------- Merge .dat files across subsessions ------------------------------------------%
%Merge .dat files from individual subsessions
%Create FileBase.cat.evt with start/stop timestamps (ms) of individual subsessions
cd(FileBase)
FileMergedOut = sprintf('%s.dat', FileBase);
fprintf('-------------------------------------------------------------------------------------------------------\n')
fprintf('            mergedat: merging subsession .dat files ---> %s \n', FileMergedOut)
fprintf('-------------------------------------------------------------------------------------------------------\n')

%Skip if the file exists
if ~exist(FileMergedOut, 'file')
    tic
    mergedat(Files2Merge, FileMergedOut)
    toc
else
    fprintf('File %s already exists. Merging was skipped. \n', FileMergedOut)
end


cd('..')
diary off
return









%---------------------------------- NOT USED ------------------------------------------------------------%
%----------------------------- Create .xml file ----------------------------------------------------------%
%Always delete the already exsting .xml file to avoid any mismath 
%between the .xml and the merged .dat files
% if exist([FileBase '.xml'], 'file')
%     delete([FileBase '.xml'])
% end
    
if ~exist([FileBase '.xml'], 'file')
    
    %Load from one abitrary .dat.info file: SamplingRate, nChan
    FileSRIn = [Files2Merge{1}, '.info'];
    if exist(FileSRIn, 'file')
        out = dlmread(FileSRIn, '\t');
        nChan = out(1);
        SamplingRate = out(2);
        clear out
    else
        error(sprintf('File %s not found!', FileSRIn))
    end
    
    
    %Load LFP sampling rate (Hz) from the parameter file
    lfpSamplingRate =  LoadMyPar(['../' filebase '/' ProcessCfgFile], 'lfpSamplingRate');
    
    %Create an .xml file with the session settings
    fprintf('Creating %s.xml ...', FileBase)
    neuroscope_xml_creator([FileBase], SamplingRate, nChan, lfpSamplingRate);
    fprintf('DONE\n')
    
end %if ~exist([FileBase '.xml'], 'file')

cd('..')










%TO DO HERE: big files convert using the chunk-approach!!!!!!
%Check whether the overall size of all the subsessions exceeds the maximum size
% if any(SubSesSize_Gb > MaxDatFileSize_Gb)
%     error(   sprintf('Size of .dat file in one or more subsessions will exceed MaxDatFileSize_Gb = %d Gb. Exclude the too large subsessions from merging.  \n', MaxDatFileSize_Gb )  )
% elseif AllSubSesSize_Gb > MaxDatFileSize_Gb
%     error(   sprintf('Size of merged .dat file will exceed MaxDatFileSize_Gb = %d Gb. Exclude the too large subsessions from merging.  \n', MaxDatFileSize_Gb )  )
% end


% %----------------- TO DO Merge multiple files from same channel if needed ----------------------------------------%
% %Check whether there are multiple files for each channel
% if length(Chan) ~= length(unique(Chan))%     
%     %Create a list of unique file names with excluded postfix
%     clear List_FileNames
%     for k=1:length(FileNames)
%         out1 = regexp(FileNames{k}, '-', 'split');
%         if length(out1)==8
%             %file name without postfix
%             out2 = regexp(FileNames{k}, '.continuous', 'split');
%         elseif length(out1)==9
%             %file name with postfix
%             out2 = regexp(FileNames{k}, '-[0-9].continuous', 'split');
%         end
%         List_FileNames{k,1} = out2{1};
%         clear out1 out2
%     end
%     List_FileNames = unique(List_FileNames);   
%     %Loop across unique file names without postfix
%     for k=1:length(List_FileNames)
%         out = dir([List_FileNames{k} '*.continuous']);
%         files2merge = arrayfun(@(x)x.name, out,'UniformOutput',0);
%         ChanFileOut = [List_FileNames{k} '-merged.continuous']% 
%         FileList = files2merge;
%         FileOut = ChanFileOut;
%         %  merge_chanfiles(files2merge, ChanFileOut)%         
%         clear out
%     end %Loop across unique file names without postfix
%     
% end  %if length(Chan) ~= length(unique(Chan)) 
% %-------------------- end ----------------------------
    
    



