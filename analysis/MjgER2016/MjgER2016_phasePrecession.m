function [P,phzStats] = MjgER2016_phasePrecession(Trial,drz,phz,spk,unit,state)


if isempty(phz),
    phz = phase(resample(Trial.load('lfp',[Trial.ephys.electrodes(:).thetaRef]),xyz),[6,12]);
end

drz = compute_drz(Trial,units);


% use unit 33 from Ed10-20140817.cof.gnd
distThresh = 250;
mResults   = 3;
fitRange   = -pi:0.01:pi;
P          = nan([numel(units),1,mResults,2]);
phzStats   = nan([numel(units),1,mResults,2]);


% GET res
res = spk(unit);
res(res>xyz.size(1))=[];


ares=res;
ares(pfdist(ares,unit==units)< distThresh) = [];
res (pfdist(res, unit==units)>=distThresh) = [];            

drzspk = drz(res,unit==units);
phzspk = phz(res,spk.map(:,unit==spk.map(:,1),2));
gind = ~isnan(drzspk)&~isnan(phzspk);

if numel(gind)<20,
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
% P:  Regression Parm: [Slope         ,Offset                          ]        
P(unit==units,s,lmid,:) = [fitRange(lmi(lmid))',atan2(sinPart(lmi(lmid)),cosPart(lmi(lmid)))'];

% Collect residuals of the theta model for each state
phi = nan([numel(lin),mResults]);
for r = lmid',
    [phi(:,r)] = polyval([2*pi;1].*sq(P(unit==units,1,r,:)),lin);
end
residuals = bsxfun(@minus,circ,phi);
residuals(residuals>pi) = residuals(residuals>pi)-2*pi;
residuals(residuals<-pi) = residuals(residuals<-pi)+2*pi;

phzStats(unit==units,s,:,:) = [circ_mean(residuals);circ_std(residuals)]';