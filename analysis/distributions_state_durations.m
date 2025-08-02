function hax = distributions_state_durations(Trial,varargin)
%function distributions_state_durations(Trial,varargin)
% 
% plots and save the log10 distribution of state durations
% 
%  varargin:
%    
%    stcMode: string, name of the state colletion which will be
%                     loaded 

defArgs = {...
    ...stcMode
       'hand_labeled_rev1',                         ...
    ...
    ...states
       {'walk','rear','turn','pause','groom','sit'},...
    ...
    ...bins
       linspace(-1,2,30),                           ...
    ...
    ...figDir
       '/storage/gravio/figures/'                   ...
};

[stcMode,states,bins,figDir] = DefaultArgs(varargin,defArgs);

% If Trial is not a MTASession try loading it.
if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end

Stc = Trial.load('stc',stcMode);

for s = Stc(states{:}),
    % unwrap state
    s = s{1}; 
    
    % plot histogram of state durations
    hfig = figure(843848);clf
    dur = diff(s.data,1,2)./s.sampleRate;
    if isempty(dur),
        hax = axes;
    else
        hax = bar(bins,histc(log10(dur),bins),'histc');
    end
    ylabel('count');
    xlabel('log10(seconds)');
    title(['Duration - ' s.label]);

    reportfig(figDir,hfig,...
              'durations',...
              'distributions',...
              false,...
              [Trial.filebase '.' s.label],...
              ['The duration distribution of ' s.label],...
              200,false,'png',4,4,[],strcmp(s.label,states{1}));
end

    