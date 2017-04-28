function Stc = get_state_transitions(Stc,Trial,varargin)
% function get_state_transitions(Stc,varargin)
% 
% Return TimePoints or TimePeriods of transitions from the
% states{1} to states{2}.
%
% INPUTS 
%     Stc:     MTAStateCollection,
%     states:  cellstr,
%     transitionWindow: numeric, 
%     ReferencData: MTAData,
%     outputType:   string, [ TimePeriods | TimePoints ], 
% 

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('states',           {{'turn','walk'}},...
                 'transitionWindow', 0.2,...
                 'ReferenceData',    [],...%Trial.load('xyz'),...
                 'outputType',       'TimePeriods'...
);

[states, transitionWindow, ReferenceData, outputType] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

smat = stc2mat(Stc,ReferenceData,states);
indexShift = round(transitionWindow * ReferenceData.sampleRate / 2);

stsTransitionMat = logical(reshape(permute(cat(3,...
                                       circshift(smat,indexShift),...
                                       circshift(smat,-indexShift)),...
                                   [1,2,3]),...
                           size(smat,1),[]));

stsTargetMat     = repmat([1,0,0,1],size(smat,1),1);

stsTransitionPeriods = ThreshCross(all(stsTransitionMat==stsTargetMat,2),0.5,0);

switch outputType,
  case 'TimePeriods'
    if nargout==0,
        Stc.addState(Trial.spath,...                         % system path to file
                     Trial.filebase,...                % prefix for filename generation
                     stsTransitionPeriods,...          % Data periods
                     ReferenceData.sampleRate,...      % Data SampleRate
                     Trial.sync.copy,...               % Trial Synchronization Periods
                     Trial.sync.data(1),...            % Synchronization Origin
                     ['t_' Stc{states{1}}.label '2' Stc{states{2}}.label],...% label
                     'T');                             % key
    else
        Stc = stsTransitionPeriods;
    end
    
  case 'TimePoints'

    
  otherwise
    error('MTA:MTAStateCollection:get_state_transitions:InvalidOutputType',...
          'You are doing it wrong.',outputType);
end



%---------------------------------------------------------------------------------------------------

