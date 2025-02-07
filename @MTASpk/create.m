function Spk = create(Spk,Session,varargin)
% Spk = create(Spk,Session,varargin)
% Function used to load spiking data from {clu,res,fet,spk} file
% and synchronize with the Session
%   
%   varargin:
%     [sampleRate,states,units,mode,loadField] 
%
%     sampleRate: double, the sampleRate in Hertz
%
%     states:     string, state label or expression to select for a
%                         subset of spikes
%
%     units:       array, Spike Cluster Ids, ([] returns all ids).
%
%     mode:       string, {'','deburst'}, only keep 1st spk of each
%                         burst.
%
%     loadField: cellarry, NOT IMPLEMENTED
%
% See also MTAStateCollection for more information on selecting
% states
%
%% Load and resample Res

% DEFARGS ----------------------------------------------------------------------
defargs = struct('sampleRate',  1,                                           ...
                     'states', [],                                           ...
                      'units', [],                                           ...
                       'mode', '',                                           ...
                  'loadField', {{}}                                          ...
);
[sampleRate, states, units, mode, loadField] =                               ...
    DefaultArgs(varargin,defargs,'--struct');
%-------------------------------------------------------------------------------


% MAIN -------------------------------------------------------------------------
filebase = fullfile(Session.spath, Session.name);

[Res, Clu, Map] = LoadCluRes( filebase );

% SELECT specific units
if isempty(units)
    cind = true( [numel(Res),1] );
elseif isnumeric(units)
    cind = find( ismember( Clu, units ));
elseif ischar(units)
    units = get_unit_set( Spk, Session, units);
    cind = find( ismember( Clu, units ));
end            
Res = Res(cind);
Clu = Clu(cind);

switch mode
    
% REMOVE burst tail spikes
  case 'deburst'
    nRes = [];
    nClu = [];
    thresh = round(0.008*Session.sampleRate);
    for u = unique(Clu)'
        try                            
            tRes   = Res( Clu==u );
            burstResInd = SplitIntoBursts( tRes, thresh);
            nRes   = [nRes; tRes(burstResInd)];
            nClu   = [nClu; u.*ones(numel(burstResInd),1)];
        end
    end
    Res = nRes;
    Clu = nClu;

  case 'blge2'
    nRes = [];
    nClu = [];
    thresh = round(0.008*Session.sampleRate);
    blThresh = 2;    
    for u = unique(Clu)'
        try
            tRes   = Res(Clu==u);
            [burstResInd,burstLength] = SplitIntoBursts(tRes,thresh);
            nRes   = [nRes; tRes(burstResInd(burstLength>=blThresh))];
            nClu   = [nClu; u.*ones([sum(burstLength>=blThresh),1])];
        end        
    end
    Res = nRes;
    Clu = nClu;    

  case 'blge3'
    nRes = [];
    nClu = [];
    thresh = round(0.008*Session.sampleRate);
    blThresh = 3;
    for u = unique(Clu)'
        try
            tRes   = Res(Clu==u);
            [burstResInd,burstLength] = SplitIntoBursts(tRes,thresh);
            nRes   = [nRes; tRes(burstResInd(burstLength>=blThresh))];
            nClu   = [nClu; u.*ones([sum(burstLength>=blThresh),1])];
        end        
    end
    Res = nRes;
    Clu = nClu;

  case 'first_spike_theta'
    % get first spike of each theta cycle

  otherwise % default
    % NOTHING
end

% RESAMPLE to target sampleRate
Res = Res ./ Session.sampleRate * sampleRate;
if sampleRate ~= 1,  Res = ceil(Res);  end

% SELECT generic periods if given
if ~isempty( Spk.per )
    newRes = [];
    newClu = [];
    for cluInd = unique(Clu)'
        sper = copy(Spk.per);
        sper.data = sper.data(Spk.perInd(cluInd,:),:) ...
            ./sper.sampleRate                         ...
            .*sampleRate;
        
        if isempty(sper.data),  continue;  end
        
        [pRes, psind] = SelectPeriods( Res, sper.data, 'd', 1, 0);
        newRes = cat(1, newRes, pRes);
        newClu = cat(1, newClu, cluInd*ones(size(pRes)));
    end
    [Res,sind] = sort(newRes);
    Clu        = newClu(sind);
end

% FIT to synchronization periods
syncPeriods = Session.sync([1,end]);
if sampleRate ~= 1
    syncPeriods =                                     ...
        ceil(syncPeriods                              ...
         * sampleRate                                 ...
         + 1 * double( sampleRate~=1 ));

end

% $$$ for per = 1:size( syncPeriods, 1)
% $$$     
% $$$ end

[Res, ind] = SelectPeriods( Res, syncPeriods, 'd', 1, 1 );
Clu = Clu(ind);


% SELECT state periods if given
if ~isempty( states );
    if ischar( states ),
        states = [Session.stc{ states, sampleRate }.data];
        [Res, sind] = SelectPeriods( Res, states, 'd', 1, 0);
    else
        sst = states.copy;
        sst.resample(sampleRate);
        [Res, sind] = SelectPeriods(Res, sst.data, 'd', 1, 0);                   
    end
    Clu = Clu(sind);    
end

Spk.clu = Clu;
Spk.res = Res;
Spk.map = Map;
Spk.sampleRate = sampleRate;

Spk.update_hash(sampleRate,mode);

% END MAIN -----------------------------------------------------------------------------------------