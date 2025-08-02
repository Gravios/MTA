function Trial = label_theta_int(Trial,varargin)
%function Trial = labelTheta(Trial,varargin)
%[Stc,thetaChan,overwrite] = DefaultArgs(varargin,{[],1,false});
% 
% RATIONAL : 
% The CA1 theta rhythm: cooridination of multiple cyclical processes.
% The EC3 input into the lacosum molecular and the CA3 input into the radiatum 
% Multiple modes of interneuron modulation
% MSDBB, EC3-FFI, CA3-FFI, CA1-FFI
% which should be reflected in the population of neurons which are
% continuous active. The theta cells, now known to be PV interneurons, are modulated at the theta
% frequency and tend to fire at a specific theta phase. 
%
% theta -> cyclic -> interneuron phase preference -> interneuron sequences -> robust phase representation???
% 
% Theta phase estimate: LFP vs INT
% 
% time series -> chuncks of time
% order
% order -> pattern
% sum of spikes given tau? definition of tau?
%
% local ccg



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('Stc',                                         [],                              ...
                 'thetaChan',                                   1,                               ...
                 'overwrite',                                   false,                           ...
);
[Stc,thetaChan,overwrite] = DefaultArgs(varargin,defargs,'--struct');

if isempty(Stc),
    Stc = Trial.stc.copy;
end
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

sempty = isempty(Stc.gsi('i'));
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('i')) = [];
    end

% LOAD lfp 
% COMPTUTE theta phase
% GET theta periods as defined by theta/delta ratio 
% GET interneurons
% LOAD spk object contanining interneurons
    lfp = Trial.load('lfp',Sessions.thetaRefGeneral);
    phz = lfp.phase([6,12]);
    tper = Trial.stc{'t',1250};
    intclu = select_units(Trial,'int','get');
    spk = Trial.spk.create(Trial,1250,[],intclu);

% COMPUTE mean theta phase and 
    mphz = [];
    rphz = [];
    for u = 1:numel(intclu),
        resu = spk(intclu(u));
        resu = resu(WithinRanges(resu,tper.data));
        mphz(u) = circ_mean(phz(resu));
        rphz(u) = circ_r(phz(resu));    
    end
    
    ufr = Trial.load('ufr',lfp,spk,intclu(rphz>0.1),0.03);
    sufr = ufr.segs(1:120:size(ufr,1),120,0);
    sufr = reshape(permute(sufr,[2,1,3]),[],size(sufr,1)*size(sufr,3));
    
    [ics,a,w] = fastica(sufr'*sufr,'numOfIC',20);
    intcmp = [];
    for c = 1:10
        u = 0;    
        while true,
            st = u*1e5+1;
            ed = (u*1e5+1e5);
            if ed>size(ufr,1),
                ed = size(ufr,1);
            end
            sufr = ufr.segs(st:ed,120,0);
            sufr = reshape(permute(sufr,[2,1,3]),[],size(sufr,1)*size(sufr,3));
            intcmp(st:ed,c) = sufr*a(:,c)/1e6*10;
            if ed==size(ufr,1),
                disp(['* Complete ',num2str(c),' *']);
                break;
            end
            u = u+1;
        end
    end
    
    ys = [];
    for c = 1:10;
        [ys(:,:,c),fs,ts] = mtchglong(WhitenSignal(intcmp(:,c)),2^11,1250,2^10,2^10*0.875,[],[],[],[0.5,20]);
    end
    
    thfout = (2 < fs & fs < 5);% | (12 < fs & fs < 14);
    thfin = 6 < fs & fs < 11;
    tdRatio = log(mean(sq(ys(:,thfin,c)),2)); %-log(mean(sq(ys(:,thfout,c)),2));
    %figure,plot(tdRatio)
    %figure,plot(log(mean(sq(ys(:,thfin,c)),2)))
    hmmStates = zeros(size(tdRatio));
    nind = ~isnan(tdRatio) & ~isinf(tdRatio) & tdRatio~=0;
    nStates = 2;
    [hmmStates(nind), thhmm, thdec] = gausshmm(tdRatio(nind),nStates,1,0);
    data = ThreshCross(hmmStates,1.5,1);
    data = round(data/diff(ts(1:2)).*lfp.sampleRate);
end

Stc.states(Stc.gsi('i')) = [];

sync = Trial.lfp.sync.copy;
%lsync = sync.sync.copy;
%lsync.resample(Trial.lfp.sampleRate);
%data = IntersectRanges(lsync.data,data)-lsync.data(1)+1;
Stc.addState(Trial.spath,...
             Trial.filebase,...
             data,...
             Trial.lfp.sampleRate,...
             Trial.sync.copy,...
             Trial.sync.data(1),...
             'thetaint','i');

Stc{'i'}.save(1);

Stc.save(1);
Trial.stc = Stc;
Trial.save;


% END MAIN -----------------------------------------------------------------------------------------