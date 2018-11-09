function [P,phzStats,Rmax,rho] = MjgER2016_phasePrecession(Trial,varargin)
%function [P,phzStats,Rmax] = MjgER2016_phasePrecession(Trial,drz,ddz,phz,spk,units)
%
% Compute phase precession coefficients between distance restricted drz and lfp phase
% 
% Varargin:
%
%    drz:        directional rate zones
%    ddz:        directional distance zones
%    phz:        Local field potential phase
%    spk:        MTASpk object which holds spike events
%    units:      list of units for computation
%    distThresh: limit phase precession analysis to specified radius around placefield center
%    mResults:   number of best fits to be returned
%
% Output:
%
%    P - NumericMatrix[(number of units), 1, : 
%    phzStats - 
%    Rmax -
%    drzHCnt - 


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('drz',                            [],                                           ...
                 'ddz',                            [],                                           ...
                 'phz',                            [],                                           ...
                 'spk',                            [],                                           ...
                 'units',                          [],                                           ...
                 'distThresh',                     250,                                          ...
                 'mResults',                       1,                                            ...
                 'fitRange',                       -pi:0.001:pi);


[drz,ddz,phz,spk,units,distThresh,mResults,fitRange]  = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

if isempty(phz),
    error('MTA_ANALYSIS_ERROR','MTA:analysis:MjgER2016:MjgER2016_phasePrecession:PhaseAbsent');
    %phz = phase(resample(Trial.load('lfp',[Trial.ephys.electrodes(:).CA1pyrThetaRef]),drz),[6,12]);
end
if isempty(drz),
    drz = compute_drz(Trial,units);
end



% use unit 33 from Ed10-20140817.cof.gnd

P          = nan([numel(units),2,mResults]);
phzStats   = nan([numel(units),2,mResults]);
s          = 1;
Rmax       = nan([numel(units),1,mResults]);
rho        = nan([numel(units),1,mResults]);

%drzHCnt    = zeros([numel(units),10]);

% GET res

for uind = 1:numel(units);
% GET unit index
    unit = units(uind);
    
    if ~isempty(spk),
% GET spikes times of unit    
% REMOVE spikes outside of position acquisition periods
% REMOVE spikes outside of distance threshold
        res = spk(unit);
        res(res>size(drz,1))=[];
        res (abs(ddz(res, uind))>=distThresh) = [];            
% SKIP fitting parameters if too few spikes
        if numel(res)<10,
            continue;
        end        
% GET direction rate zone (DRZ) values at times of spikes ( see Huxter(2008) )
% GET phase values at times of spikes   
% IGNORE spikes where drz or phase are nans
        drzspk = drz(res,uind);
        phzspk = phz(res,spk.map(unit==spk.map(:,1),2));            
    else,% compute phase precession for single unit
        drzspk = drz(abs(ddz)<=distThresh);
        phzspk = phz(abs(ddz)<=distThresh);
    end
        
    gind = ~isnan(drzspk)&~isnan(phzspk);

% SKIP fitting parameters if too few spikes
    if sum(gind)<10,
        continue;
    end    
    

    %drzHCnt(uind,:) = histcounts(drzspk,linspace(-1,1,11));


    lin = drzspk(gind); 
    circ = phzspk(gind);
    
    cosPart = sum(cos(bsxfun(@minus,circ,2*pi*bsxfun(@times,fitRange,lin))),1);
    sinPart = sum(sin(bsxfun(@minus,circ,2*pi*bsxfun(@times,fitRange,lin))),1);
    
    R = sqrt((cosPart./length(circ)).^2+...
             (sinPart./length(circ)).^2 );

    [lmi,~] = LocalMinima(-R',0,0,mResults);

    lmid = find(~isnan(lmi));
    % Rmax fit quality 
    Rmax(uind,1,lmid) = R(lmi(lmid));
    % P:  Regression Parm: [Slope,Offset]
    P(uind,:,lmid) = permute([2*pi*fitRange(lmi(lmid))',atan2(sinPart(lmi(lmid)),cosPart(lmi(lmid)))'],[3,2,1]);

    % Collect residuals of the theta model for each state
    phi = nan([numel(lin),mResults]);
    for r = lmid',
        [phi(:,r)] = polyval([2*pi;1].*sq(P(uind,:,r))',lin);
    end
    residuals = bsxfun(@minus,circ,phi);
    residuals(residuals>pi) = residuals(residuals>pi)-2*pi;
    residuals(residuals<-pi) = residuals(residuals<-pi)+2*pi;

    phzStats(uind,:,:) = [circ_mean(residuals);circ_std(residuals)];

% COMPUTE circular-linear correlation coefficient
    linC = mod(abs(P(uind,1,1))*lin,2*pi);
    circMean = atan2(sum(sin(circ)),sum(cos(circ)));
    linCMean = atan2(sum(sin(linC)),sum(cos(linC)));
    rho(uind,1,lmid) = sum(sin(circ-circMean).*sin(linC-linCMean))...
               ./sqrt(sum(sin(circ-circMean).^2).*sum(sin(linC-linCMean).^2));
    
end

% $$$ figure,
% $$$ hold('on');
% $$$ plot([lin,lin],[circ,circ+2*pi],'b.')
% $$$ plot([-1,1],P(uind,1,1)*[-1,1]+P(uind,2,1),'-m','LineWidth',1)
% $$$ plot([-1,1],P(uind,1,1)*[-1,1]+P(uind,2,1)+2*pi,'-m','LineWidth',1)
