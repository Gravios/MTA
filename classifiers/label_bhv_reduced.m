function Stc = label_bhv_reduced(Stc,Trial,varargin);


% DEFARGS -----------------------------------------------------------------------------------------
defargs = struct('stcMode',                     'msnn_ppsvd',                                  ...
                 'tag',                         '_raux'                                        ...
);
[stcMode,tag] = DefaultArgs(varargin,defargs,'--struct');
% -------------------------------------------------------------------------------------------------

if isempty(Stc),
    Stc = Trial.load('stc',stcMode);
end

Stc = reduce_stc_to_loc(Stc);

% LABEL high and low versions of locomotion and pause
Stc = label_bhv_reduced_aux(Stc,Trial);

% SET modified name for the state collection mode
Stc.updateMode([Stc.mode,tag]);

Stc.save(1);