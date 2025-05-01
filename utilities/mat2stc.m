function [Stc] = mat2stc(mstc,Stc,Data,Trial,varargin)
%function [Stc] = stc2mat(mstc,Stc,Data)
% assumes no heirarchical relationships between states
% each is mutually exclusive from all others

    
% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(...
    'labels',             {Stc.list_state_attrib('label')},...
    'keys',               {Stc.list_state_attrib('key')}...
);
[labels,keys] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

for i = 1:numel(labels),
    sts = Stc.gsi(labels{i});
    if isempty(sts),
        Stc.addState(Trial.spath,...
                     Trial.filebase,...
                     mstc(:,i)>0,Data.sampleRate,...
                     Trial.sync.copy,...
                     Trial.sync.data(1),...
                     labels{i},...
                     keys{i},...
                     'type','TimeSeries');
        cast    (Stc.states{sts},'TimePeriods');
    else        
        resample(Stc.states{sts}, Data);
        cast    (Stc.states{sts}, 'TimeSeries');
        Stc.states{sts}.data = mstc(:,i)>0;
        cast    (Stc.states{sts}, 'TimePeriods');
    end    

    clean   (Stc.states{sts});
    resample(Stc.states{sts},Data);        
end

% END MAIN -----------------------------------------------------------------------------------------