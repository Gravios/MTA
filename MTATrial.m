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
            Trial.filebase = [Trial.name '.' Trial.maze.name '.' Trial.trialName];
            Trial.stc.updateFilename(Trial.filebase);
            
            
            if exist(fullfile(Trial.spath, [Trial.filebase '.trl.mat']),'file')&&~overwrite
                ds = load(fullfile(Trial.spath, [Trial.filebase '.trl.mat']));
                %Trial.xyzPeriods = ds.xyzPeriods;
                new_xyzPeriods = ds.xyzPeriods;
                if isfield(ds,'bhvmode'),
                    if ~isempty(ds.bhvmode)&&exist(fullfile(Trial.spath, [Trial.filebase '.stc.' ds.bhvmode '.mat']),'file')
                        Trial.stc.updateFilename([Trial.filebase '.stc.' ds.bhvmode '.mat']);
                        Trial.stc.load;
                        Trial.stc.updateSync(Trial.sync);
                    end
                end
                if strcmp(mode,'minimal'),
                    Trial.xyz = MTADxyz(Trial.spath,Trial.name,[],[]);
                    return,
                end
            end
            
            Trial.trackingMarker = Session.trackingMarker;
            Trial.sync.resync(Trial.xyz,new_xyzPeriods);
            props = properties(Trial);
            for p = 1:numel(props),
                prop = Trial.(props{p});
                if isa(prop,'MTAData'),
                    if prop.isempty||strcmp(props{p},'xyz'),
                        continue,
                    else
                        Trial.sync.resync(prop);
                    end
                elseif isa(prop,'MTAStateCollection')
                    if ~prop.isempty
                        for s = numel(prop(:)),
                            Trial.sync.resync(prop{s});
                        end
                    end
                end
            end
            
            if numel(preLoadedFields)>0;
                for f = 1:numel(preLoadedFields),
                    field = preLoadedFields{f};
                    if iscell(field)
                        Trial.(field{1}).load(field(2:end))
                        if ~Trial.(field{1}).isempty,
                            Trial.sync.resync(Trial.(field{1}));
                        end
                    else
                        Trial.(field).load;
                        if ~Trial.(field).isempty,
                            Trial.sync.resync(Trial.(field));
                        end
                    end
                end
            end
        end
        

        function save(Trial)
            bhvmode = [];
            if ~isempty(Trial.stc)
                bhvmode = Trial.stc.mode;
            end
            trialName = Trial.trialName;
            xyzPeriods = Trial.sync.periods(Trial.xyz.sampleRate);
            save(fullfile(Trial.spath,[Trial.filebase '.trl.mat']),'trialName','xyzPeriods','bhvmode');            
        end

    end
end
