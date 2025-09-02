function Spk = create(Spk,Session)
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
% TODO : change name to generate
%

% >>> DEFARGS >>> -------------------------------------------------------------
inp = inputParser();

inp.addRequired ( 'Spk'                );
inp.addRequired ( 'Session'            );

inp.addParameter(  'sampleRate', 1     );
inp.addParameter(      'states', []    );
inp.addParameter(       'units', 'all' );
inp.addParameter(        'mode', ''    );
inp.addParameter(   'incld_fet', false );
inp.addParameter(   'incld_swf', false );
inp.parse(varargin{:});
[Spk, Session, sampleRate, states, units, mode, incld_fet, incld_swf] = deal ...
    (                                                                        ...
        inp.Results.Spk,                                                    ...
        inp.Results.Session,                                                ...
        inp.Results.sampleRate,                                             ...
        inp.Results.states,                                                 ...
        inp.Results.units,                                                  ...
        inp.Results.mode,                                                   ...
        inp.Results.incld_fet,                                              ...
        inp.Results.incld_swf                                               ...
    );
% <<< DEFARGS <<< -------------------------------------------------------------


% >>> MAIN >>> ----------------------------------------------------------------

Spk.load_clu_res();

% >>> SELECT specific units >>> -----------------------------------------------
if isempty(units)
    cind = true( [numel(Spk.res),1] );
elseif isnumeric(units)
    cind = find( ismember( Spk.clu, units ));
elseif ischar(units)
    units = get_unit_set( Spk, Session, units);
    cind = find( ismember( Spk.clu, units ));
end            
Spk.res = Spk.res(cind);
Spk.clu = Spk.clu(cind);
if incld_fet, Spk.fet = Spk.fet(cind,:); end
if incld_spk, Spk.spk = Spk.spk(cind,:,:);end    
% >>> SELECT specific units >>> -----------------------------------------------

% >>> FILTER mode >>> ---------------------------------------------------------
switch mode
  case 'deburst'
  % >>> REMOVE burst tail spikes >>> ------------------------------------------
    nRes = [];
    nClu = [];
    nFet = [];
    nSpk = [];
    thresh = round(0.008*Session.sampleRate);
    for u = unique(Spk.clu)'
        try                            
            tRes   = Spk.res( Spk.Clu==u );
            burstResInd = SplitIntoBursts( tRes, thresh);
            nRes   = [nRes; tRes(burstResInd)];
            nClu   = [nClu; u.*ones(numel(burstResInd),1)];
            if incld_fet,  nFet = [nFet; tFet(burstResInd,:)];   end
            if incld_spk,  nSpk = [nSpk; tSpk(burstResInd,:,:)]; end
        end
    end
    Spk.res = nRes;
    Spk.clu = nClu;
    if incld_fet,  Spk.fet = nFet; end
    if incld_spk,  Spk.spk = nSpk; end
  % <<< REMOVE burst tail spikes <<< ------------------------------------------  
  % >>> other modes >>> -------------------------------------------------------
% $$$   case 'blge2'
% $$$     nRes = [];
% $$$     nClu = [];
% $$$     thresh = round(0.008*Session.sampleRate);
% $$$     blThresh = 2;    
% $$$     for u = unique(Clu)'
% $$$         try
% $$$             tRes   = Spk.res(Spk.clu==u);
% $$$             [burstResInd,burstLength] = SplitIntoBursts(tRes,thresh);
% $$$             nRes   = [nRes; tRes(burstResInd(burstLength>=blThresh))];
% $$$             nClu   = [nClu; u.*ones([sum(burstLength>=blThresh),1])];
% $$$         end        
% $$$     end
% $$$     Spk.res = nRes;
% $$$     Spk.clu = nClu;    
% $$$ 
% $$$   case 'blge3'
% $$$     nRes = [];
% $$$     nClu = [];
% $$$     thresh = round(0.008*Session.sampleRate);
% $$$     blThresh = 3;
% $$$     for u = unique(Clu)'
% $$$         try
% $$$             tRes   = Spk.res(Spk.clu==u);
% $$$             [burstResInd,burstLength] = SplitIntoBursts(tRes,thresh);
% $$$             nRes   = [nRes; tRes(burstResInd(burstLength>=blThresh))];
% $$$             nClu   = [nClu; u.*ones([sum(burstLength>=blThresh),1])];
% $$$         end        
% $$$     end
% $$$     Spk.res = nRes;
% $$$     Spk.clu = nClu;
  % <<< other modes <<< -------------------------------------------------------
  case 'first_spike_theta'
    % get first spike of each theta cycle
  otherwise % default
    % NOTHING for now
end
% <<< FILTER mode <<< ---------------------------------------------------------

% >>> RESAMPLE to target sampleRate >>> ---------------------------------------
Spk.res = Spk.res ./ Session.sampleRate * sampleRate;
if sampleRate ~= 1,  Spk.res = ceil(Spk.res);  end
Spk.sampleRate = sampleRate;
% <<< RESAMPLE to target sampleRate <<< ---------------------------------------

% >>> FIT to Trial synchronization periods >>> --------------------------------
syncPeriods = Session.sync([1,end]);
if sampleRate ~= 1
    syncPeriods = ceil(syncPeriods * sampleRate + 1 * double( sampleRate~=1 ));
end
[Spk.res, ind] = SelectPeriods( Spk.res, syncPeriods, 'd', 1, 1 );
Spk.clu = Spk.clu(ind);
if incld_fet,  Spk.fet = Spk.fet(ind,:);   end
if incld_spk,  Spk.spk = Spk.spk(ind,:,:); end
% <<< FIT to Trial synchronization periods <<< --------------------------------

% >>> SELECT state periods if given >>> ---------------------------------------
if ~isempty( states );
    if ischar( states ),
        states = [Session.stc{ states, sampleRate }.data];
        [Spk.res, sind] = SelectPeriods( Spk.res, states, 'd', 1, 0);
    else
        sst = states.copy;
        sst.resample(sampleRate);
        [Spk.res, sind] = SelectPeriods(Spk.res, sst.data, 'd', 1, 0);
    end
    Spk.clu = Spk.clu(sind);
    if incld_fet,  Spk.fet = Spk.fet(sind,:);   end
    if incld_spk,  Spk.spk = Spk.spk(sind,:,:); end
end
% <<< SELECT state periods if given <<< ---------------------------------------

Spk.update_hash(sampleRate,mode);

% <<< MAIN <<< ----------------------------------------------------------------
