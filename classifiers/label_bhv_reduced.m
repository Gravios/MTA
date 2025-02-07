function Stc = label_bhv_reduced(Stc,Trial,varargin);


% DEFARGS -----------------------------------------------------------------------------------------
defargs = struct('tag','_raux');
[tag] = DefaultArgs(varargin,defargs,'--struct');
% -------------------------------------------------------------------------------------------------

if ischar(Stc),
    Stc = Trial.load('stc',Stc);
end

Stc = reduce_stc_to_loc(Stc, Trial);

% LABEL high and low versions of locomotion and pause
Stc = label_bhv_reduced_aux(Stc, Trial);

% SET modified name for the state collection mode
Stc.updateMode([Stc.mode,tag]);

Stc.save(1);