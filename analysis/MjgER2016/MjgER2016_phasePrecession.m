function [P,phzStats,Rmax] = MjgER2016_phasePrecession(Trial,drz,ddz,phz,spk,units)
%function [P,phzStats,Rmax] = MjgER2016_phasePrecession(Trial,drz,ddz,phz,spk,units)
%
% Compute phase precession coefficients between distance restricted drz and lfp phase
% 
% Output:
%
%    P - NumericMatrix[(number of units), 1, : 

if isempty(phz),
    error('MTA_ANALYSIS_ERROR','MTA:analysis:MjgER2016:MjgER2016_phasePrecession:PhaseAbsent');
    %phz = phase(resample(Trial.load('lfp',[Trial.ephys.electrodes(:).CA1pyrThetaRef]),drz),[6,12]);
end
if isempty(drz),
    drz = compute_drz(Trial,units);
end


% use unit 33 from Ed10-20140817.cof.gnd
distThresh = 250;
mResults   = 3;
fitRange   = -pi:0.01:pi;
P          = nan([numel(units),mResults,2]);
phzStats   = nan([numel(units),mResults,2]);
s          = 1;
Rmax       = nan([numel(units),mResults]);

% GET res

for unit = units;
% GET unit index
    uind = find(unit==units);
    
% GET spikes times of unit    
% REMOVE spikes outside of position acquisition periods
% REMOVE spikes outside of distance threshold
    res = spk(unit);
    res(res>size(drz,1))=[];
    ares=res;
    ares(ddz(ares,uind)< distThresh) = [];
    res (ddz(res, uind)>=distThresh) = [];            

% SKIP fitting parameters if too few spikes
    if numel(res)<10,
        continue;
    end    
    
% GET direction rate zone (DRZ) values at times of spikes ( see Huxter(2008) )
% GET phase values at times of spikes   
% IGNORE spikes where drz or phase are nans
    drzspk = drz(res,uind);
    phzspk = phz(res,spk.map(unit==spk.map(:,1),2));    
    gind = ~isnan(drzspk)&~isnan(phzspk);

% SKIP fitting parameters if too few spikes
    if sum(gind)<10,
        continue;
    end    


    lin = drzspk(gind); 
    circ = phzspk(gind);
    x=[min(lin),max(lin)];

    cosPart = sum(cos(bsxfun(@minus,circ,2*pi*bsxfun(@times,fitRange,lin))),1);
    sinPart = sum(sin(bsxfun(@minus,circ,2*pi*bsxfun(@times,fitRange,lin))),1);
    R = sqrt((cosPart./length(circ)).^2+...
             (sinPart./length(circ)).^2 );


    [lmi,lmv] = LocalMinima(-R',0,0,mResults);

    lmid = find(~isnan(lmi));
    % Rmax fit quality 
    Rmax(uind,lmid) = R(lmi(lmid));
    % P:  Regression Parm: [Slope,Offset]
    P(uind,lmid,:) = [fitRange(lmi(lmid))',atan2(sinPart(lmi(lmid)),cosPart(lmi(lmid)))'];

    % Collect residuals of the theta model for each state
    phi = nan([numel(lin),mResults]);
    for r = lmid',
        [phi(:,r)] = polyval([2*pi;1].*sq(P(uind,r,:)),lin);
    end
    residuals = bsxfun(@minus,circ,phi);
    residuals(residuals>pi) = residuals(residuals>pi)-2*pi;
    residuals(residuals<-pi) = residuals(residuals<-pi)+2*pi;

    phzStats(uind,:,:) = [circ_mean(residuals);circ_std(residuals)]';

end
