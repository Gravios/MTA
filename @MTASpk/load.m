function Spk = load(Spk, Session, varargin)
% Spk = load(Spk,Session,varargin)
% Function used to load spiking data from {clu,res,fet,spk} file
% and synchronize with the Session
%   
%   varargin:
%     [Spk,sampleRate,states,units,mode,loadField] 
%
%     sampleRate: double, the sampleRate in Hertz
%
%     states:     string, state label or expression to select for a
%                         subset of spikes
%
%     units:       array, Cluster identities to select for a subset
%                         of spikes. The default [] will return all
%                         Clusters.
%
%     mode:       string, {'','deburst'}, only keep 1st spk of each
%                         burst.
%
%     incld_fet:   Logical, {false}, load the spike PCA Features
%
%     incld_swf:   Logical, {false}, load the spike wave forms
%
%     denoise:     Logical, {true}, remove artifact and noise clusters
%
% See MTAStateCollection for more information on selecting states
%

% >>> DEFARGS >>> -------------------------------------------------------------
defargs = struct( 'sampleRate' , 1,                                         ...
                      'states' , [],                                        ...
                       'units' , 'all',                                     ...
                        'mode' , '',                                        ...
                   'incld_fet' , false,                                     ...
                   'incld_swf' , false,                                     ...
                     'denoise' , true                                       ...
                );
[sampleRate, states, units, mode, incld_fet, incld_swf, denoise] =          ...
    DefaultArgs(varargin, defargs, '--struct');
% <<< DEFARGS <<< -------------------------------------------------------------

% >>> MAIN >>> ----------------------------------------------------------------

load_ndm_spikes(Spk, units, incld_fet, incld_swf, denoise);

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
            if incld_swf,  nSpk = [nSpk; tSpk(burstResInd,:,:)]; end
        end
    end
    Spk.res = nRes;
    Spk.clu = nClu;
    if incld_fet,  Spk.fet = nFet; end
    if incld_swf,  Spk.spk = nSpk; end
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
    syncPeriods = ceil(syncPeriods .* sampleRate + 1 .* double( sampleRate~=1 ));
end
[Spk.res, ind] = SelectPeriods( Spk.res, syncPeriods, 'd', 1, 1 );
Spk.clu = Spk.clu(ind);
if incld_fet,  Spk.fet = Spk.fet(ind,:);   end
if incld_swf,  Spk.spk = Spk.spk(ind,:,:); end

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
    if incld_swf,  Spk.spk = Spk.spk(sind,:,:); end
end
% <<< SELECT state periods if given <<< ---------------------------------------

Spk.update_hash( sampleRate, mode);

% <<< MAIN <<< ----------------------------------------------------------------

