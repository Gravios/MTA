function Data = load(Data,Session,varargin)
%load(Data,Session,varargin)
% feature    
%[channels,gselect,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},[]});
%
% fiberFiles = dir(['ZY0519-20190820.fbr.*'])
% fiberFiles = dir(['ZY0519-20190820.fbr.hpc.*'])
% fiberFiles = dir(['ZY0519-20190820.fbr.hpcCA1.NULL'])
% fiberFiles = dir(['ZY0519-20190820.fbr.hpcCA1.Dopa'])
% fiberFiles = dir(['ZY0519-20190820.fbr.hpcCA1.ACh'])
% fiberFiles = dir(['ZY0519-20190820.fbr.hpcCA1.N'])
% % LOAD all signals for hpcCA1
% ACh = Session.load('fbr','hpcCA1');
% 
% % LOAD all ACh signal for all locations
% ACh = Session.load('fbr',[],'ACh');

    
    
%%%<<< DESCRIPTION
% LOAD cell array of structures containing the fiber data and features
%     NOTE : fiber data is synced to lfp leading edge is at the lfp origin
% 
% SET origin to 1/(sampleRate/2)
%%%>>>

%%%<<< TESTING VARS
% $$$ location = 'hpc';
% $$$ signalType = 'SENR';
% $$$ field = 'signal';
% $$$ Session.name = 'ZY0519-20190820';
% $$$ Session.spath = '/storage2/ziyan/data/processed/spect/ZY0519/ZY0519-20190820/';
%%%>>>    
    
    
% DEFARGS ------------------------------------------------------------------------------------------    

[location,signalType,field,periods] = DefaultArgs(varargin,{'','','signal',[]});
%---------------------------------------------------------------------------------------------------



% LOAD data from *.fbr 
%filename = 'ZY0519-20190820.fbr.hpc.SENR.mat';
filename = [Session.name,'.fbr.',location,'.',signalType,'.mat'];
fbr = load( fullfile( Session.spath(), filename));

%fbr = load( fullfile( Session.spath(), filename),...
%           'sampleRate', 'intergrationTime', field);
%
%           sigHbClean: [10995×845 double]
%                hbFit: [10995×3 double]
%              sigNorm: [10995×845 double]
%                    t: [10995×1 double]
%             spectACh: [10995×845 double]
%               signal: [10995×1 double]
%           sampleRate: 10
%     intergrationTime: 0.09
%             location: 'hpc'
%           signalType: 'SENR'


% SET synchronization periods
if isempty(periods),
    Session.sync.resample(1);
    periods = round(Session.sync([1,end]).*fbr.sampleRate);
elseif isa(periods,'MTADepoch'),
    periods.resample(Session.fbr.sampleRate);
    periods = periods.data + Session.sync(1)*fbr.sampleRate;
else
    periods = periods*fbr.sampleRate; % ADD GENERAL PERIOD SELECTION
end
periods(periods==0) = 1;
periods(periods>numel(fbr.signal)) = numel(fbr.signal);

Data.clear();

name = strjoin({location,signalType,field},'_');
Data.integrationTime = fbr.integrationTime;
Data.model = [];
Data.data  = [];
Data.path  = Session.spath;
Data.filename = filename;
Data.name  = ['fiber_photometry_',name];
Data.label = name;
Data.key   = 'f';
Data.sampleRate = fbr.sampleRate;
Data.sync  = copy(Session.lfp.sync);
Data.origin = 0;

for period = periods',
    Data.data = cat(1,Data.data,fbr.(field)(period(1):period(2),:,:));
end

