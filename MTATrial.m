classdef MTATrial < MTASession
% MTATrial(name,varargin) - Is a subclass of MTASession structured to organize 
%   subsets of a parent session to aid in the analysis of neural and spatial data.
%
%   Ext: *.trl.mat
%   Saved Trial object contain only basic information required to load subsets of
%   a full session
% 
%   name - string: Same name as the directory of the session
%   
%   -or- 
%
%   name - MTASession: The Session object from which data will be copied
%
%
%   varargin:
%     [preLoadedFields,trialName,new_xyzPeriods,overwrite,mazeName]
%    
%     preLoadedFields: cellArray, load saved data into Trial fields
%                      (e.g. {'ufr','ang'} or {'ufr',{'lfp',65:96},'CluRes'}
%
%     trialName:       string, designation of trial, the Trial representing 
%                      the full Session has the default name 'all'
%
%     new_xyzPeriods:  numericArray, new Periods defined in the refrence 
%                      frame of the MTASession objects xyzPeriods
%
%     overwrite:       boolean, flag to overwrite saved Trials
%
%     mazeName:        string, 3-4 letter name of the testing arena 
%                      (e.g. 'rof' := rectangular open-field)
%
%
%---------------------------------------------------------------------------------------------------------
%   General Loading:
%     
%     Load from saved Trial,
%     Trial = MTATrial(name,[],trialName);
%     
%     Create new Trial,
%     Trial = MTATrial(name,[],trialName,new_xyzPeriods,overwrite,mazeName);
%
%
%---------------------------------------------------------------------------------------------------------
%     examples:
%       load saved Trial,
%         Trial = MTATrial('jg05-20120309',[],'all');
%
%       load saved Trial with presaved variables: lfp and Place fields
%         Trial = MTATrial('jg05-20120309',{{'lfp',65:95},'Pfs'},'all');
%
%       load saved Trial with presaved variables: marker angles and Clustering data
%         Trial = MTATrial('jg05-20120309',{'ang','CluRes'},'all');
%
%       Create New Trial from a subset of the total session
%         Trial = MTATrial('jg05-20120309',[],'crt1',[1,10000;11000,21000],1,'rof','normal')
%
%              
%---------------------------------------------------------------------------------------------------------

    methods 
        function Trial = MTATrial(Session,varargin)
            [preLoadedFields,trialName,new_xyzPeriods,overwrite,mazeName,mode] = DefaultArgs(varargin,{{},'all',[],0,'cof','normal'});
            if ~isa(Session,'MTASession'),
                Session = MTASession(Session,{},mazeName);
            end
            if strcmp(trialName,'nil'),overwrite = 1;end
            
            Trial = Trial@MTASession(Session,{},mazeName);
            
            Trial.trialName = trialName;
            Trial.filebase = [Trial.name '.' Trial.Maze.name '.' Trial.trialName];
            Trial.Stc = {};
            if exist(fullfile(Trial.spath, [Trial.filebase '.trl.mat']),'file')&&~overwrite
                ds = load(fullfile(Trial.spath, [Trial.filebase '.trl.mat']));
                Trial.xyzPeriods = ds.xyzPeriods;
                if isfield(ds,'bhvmode'),
                    if ~isempty(ds.bhvmode)&&exist(fullfile(Trial.spath, [Trial.filebase '.Stc.' ds.bhvmode '.mat']),'file')
                        Trial.Stc.load(ds.bhvmode);
                    end
                end
                if strcmp(mode,'minimal'),
                    Trial.xyz=[];
                    Trial.xyzSegLength=[];
                    return,
                end
            elseif ~ischar(new_xyzPeriods)
                Trial.xyzPeriods = new_xyzPeriods;
            else
                switch new_xyzPeriods
                  case 'default_height_restricted'
                    height_violations = [];
                    try
                        height_violations = ThreshCross(Trial.xyz(:,1,3),300,10);
                    end
                    if ~isempty(height_violations),
                        for i = 1:size(height_violations,1),
                            for j = 1:size(Trial.xyzPeriods,1),
                                if height_violations(i,2)>Trial.xyzPeriods(j,1)&height_violations(i,2)<(Trial.xyzPeriods(j,1)+3000),
                                    Trial.xyzPeriods(j,1) = height_violations(i,2)+400;
                                elseif height_violations(i,1)<Trial.xyzPeriods(j,2)&height_violations(i,1)>(Trial.xyzPeriods(j,2)-3000),
                                    Trial.xyzPeriods(j,2) = height_violations(i,1)-400;
                                else
                                    continue
                                end
                            end
                        end
                    end
                end
            end

            % Reset other properties with new xyzPeriods            
            Trial.xyzSegLength = diff(Trial.xyzPeriods,1,2)+1;
            Trial.trackingMarker = Session.trackingMarker;
            syncShift = Trial.syncPeriods(1,1);
            Trial.syncPeriods = [];
            Trial.syncPeriods = round((Trial.xyzPeriods-1.05)/Trial.xyz.sampleRate*Trial.lfpSampleRate+syncShift);

            xyz = false(Trial.xyz.size(1),1);            
            for i = 1:size(Trial.xyzPeriods,1),
                xyz(Trial.xyzPeriods(i,1):Trial.xyzPeriods(i,2))=true;               
            end
            Trial.xyz.data(xyz==0,:,:) = 0;            
            Trial.xyz.data = Trial.xyz(Trial.xyzPeriods(1):Trial.xyzPeriods(end),:,:);

            % resync all non-empty fields
            props = properties(Trial);
            for p = 1:numel(props),
                if isa(Trial.(props{p}),'MTAData'),
                    if Trial.(props{p}).isempty||strcmp(props{p},'xyz'),
                        continue,
                    else
                        if Trial.(props{p}).sampleRate==Trial.xyz.sampleRate,
                            Trial.(props{p}).data(xyz==0,:,:,:,:) = 0;
                            Trial.(props{p}).data = Trial.(props{p})(Trial.xyzPeriods(1):Trial.xyzPeriods(end),:,:,:,:);
                        elseif Trial.(props{p}).sampleRate==Trial.lfp.sampleRate,
                            Trial.(props{p}).data = Trial.(props{p}).data(Trial.syncPeriods(1):Trial.syncPeriods(end),:,:,:,:);
                        end
                    end
                end
            end
            
            % resync the stc.states
            %for s=numel(s.stc.states)
            
            
            
            %Load fields such as ang,lpf,spk,...ect. within the new time
            %periods.
            if numel(preLoadedFields)>0;
                for f = 1:numel(preLoadedFields),
                    field = preLoadedFields{f};
                    %Why is this try here?
                    try
                        if iscell(field)
                            Trial.(field{1}).load(field(2:end))
                            if ~isempty(Trial.(field{1}).data),
                                if Trial.xyz.sampleRate==Trial.(field{1}).sampleRate
                                    Period = Trial.xyzPeriods(1,end);
                                else
                                    Period = Trial.syncPeriods(1,end);
                                end
                                Trial.(field{1}).data = Trial.(field{1}).data(Period);
                            end
                        else
                            Trial.(field).load;
                            if ~isempty(Trial.(field).data),
                                if Trial.xyz.sampleRate==Trial.(field).sampleRate
                                    Period = Trial.xyzPeriods([1,end]);
                                else
                                    Period = Trial.syncPeriods([1,end]);
                                end
                                Trial.(field).data = Trial.(field).data(Period);
                            end
                        end
                    catch
                        
                    end                    
                end
            end
        end
        

        function save(Trial)
            bhvmode = [];
            if ~isempty(Trial.Bhv)
                bhvmode = Trial.Bhv.mode;
            end
            trialName = Trial.trialName;
            xyzPeriods = Trial.xyzPeriods;
            save(fullfile(Trial.spath,[Trial.filebase '.trl.mat']),'trialName','xyzPeriods','bhvmode');            
        end

    end
end
