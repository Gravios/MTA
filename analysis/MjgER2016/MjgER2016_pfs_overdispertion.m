
%p(s|x,y)
%p(s|x,y,h,b)

global MTA_PROJECT_PATH
% LOAD analysis data
MjgER2016_load_data();
% SETUP default args for analysis functions
MjgER2016_pfs_overdispertion_args('section 1');

sampleRate = 30; % Hz
ufrBinWidth = 2; % seconds
interpParPfs = struct('bins',{{linspace(-500,500,50),    ... mm
                               linspace(-500,500,50),    ... mm
                               linspace(  -2,  0.8,50),  ... rad
                               linspace(  -0.8,2,  50)}},... rad
                      'nanMaskThreshold', 0.1,           ...
                      'methodNanMap',     'linear',      ...
                      'methodRateMap',    'linear');
% LOAD ratemap objects
pft       = cf(@(t,u)   pfs_2d_theta(t,u,'pfsArgsOverride',...
                                    struct('halfsample',false,'numIter',1)), ...
                        Trials, units);
ddz       = cf(@(t,u,p) compute_ddz(t,u,p,'sampleRate',sampleRate),          ...
                        Trials, units,pft);
% GET list of well sampled unitsn
bfs       = cf(@(t,u)   compute_bhv_ratemaps(t,u),   Trials, units);
[LR,FSr,VT,unitSubset,validDims,~,~] = compute_bhv_ratemaps_erpPCA(bfs,units);
cluSessionMapSubset = cluSessionMap(unitSubset,:);
units = cell([1,numel(Trials)]);
for unit = cluSessionMapSubset',
    units{unit(1)} = [units{unit(1)},unit(2)];
end

pfs.xy    = cf(@(t,u)   compute_ratemaps(t,u),       Trials, units);
pfs.xyhb  = cf(@(t,u)   compute_xyhb_ratemaps(t,u),  Trials, units);


% LOAD position and behavior space objects
xyz = cf(@(t)  resample(preproc_xyz(t,'trb'),sampleRate),   Trials);
      cf(@(x)  set(x,'data',clip(x.data,-499,499)),         xyz   );
fet = cf(@(t)  fet_HB_pitchB(t,sampleRate),                 Trials);
      cf(@(f)  set(f,'data',[clip(f.data(:,1),-1.99,0.799),clip(f.data(:,2),-0.799,1.99)]), fet);
% LOAD binned firing rates
ufr = cf(@(t,u,x) t.load('ufr',x,[],u,1,[],'count'), Trials,units,xyz);      
% GET indicies for each variable within the ratemap space
xind = cf(@(x) discretize(x(:,'hcom',1),linspace(-500,500,pfs.xyhb{1}.adata.binSizes(1)+1)), xyz);
yind = cf(@(x) discretize(x(:,'hcom',2),linspace(-500,500,pfs.xyhb{1}.adata.binSizes(2)+1)), xyz);
hind = cf(@(f) discretize(f(:,1),linspace(-2,0.8,pfs.xyhb{1}.adata.binSizes(3)+1)),        fet);
bind = cf(@(f) discretize(f(:,2),linspace(-0.8,2,pfs.xyhb{1}.adata.binSizes(4)+1)),        fet);
% CONVERT 
indXY   = cf(@(pf,xi,yi)       sub2ind(pf.adata.binSizes',xi,yi),       pfs.xy,   xind, yind);
indXYHB = cf(@(pf,xi,yi,hi,bi) sub2ind(pf.adata.binSizes',xi,yi,hi,bi), pfs.xyhb, xind, yind, hind, bind);

expRate.xy = cf(@(pf,unt,ixy) ...
                af(@(u) mean(pf.data.rateMap(ixy,pf.data.clu==u,:),3,'omitnan'), unt), ...
                pfs.xy, units, indXY);
expRate.xyhb = cf(@(pf,unt,ixy) ...
                  af(@(u) mean(pf.data.rateMap(ixy,pf.data.clu==u,:),3,'omitnan'), unt), ...
                  pfs.xyhb, units, indXYHB);

expRate.xy   = cf(@(e) cat(2,e{:}), expRate.xy);
expRate.xyhb = cf(@(e) cat(2,e{:}), expRate.xyhb);

expCount.xy   = cf(@(e) ...
                   sq(sum(GetSegs(e,1:sampleRate.*ufrBinWidth:size(e,1),...
                                  sampleRate.*ufrBinWidth,nan).*1./sampleRate,'omitnan')), ...
                   expRate.xy);
expCount.xyhb = cf(@(e) ...
                   sq(sum(GetSegs(e,1:sampleRate.*ufrBinWidth:size(e,1),...
                                  sampleRate.*ufrBinWidth,nan).*1./sampleRate,'omitnan')), ...
                   expRate.xyhb);
ddzMask  = cf(@(e) sq(all([any(abs(GetSegs(e,1:sampleRate.*ufrBinWidth:size(e,1),...
                                           sampleRate.*ufrBinWidth,1))<100,1);...
                    mean(abs(GetSegs(e,1:sampleRate.*ufrBinWidth:size(e,1),...
                                     sampleRate.*ufrBinWidth,1)))<200],1)), ...
              ddz);



tper = cf(@(t) [t.stc{'t-m-s',30}], Trials);
tper = cf(@(t) [t.stc{'w+r+n&t',30}], Trials);

tind = cf(@(p,e,d) WithinRanges((sampleRate.*ufrBinWidth/2+1):sampleRate.*ufrBinWidth:size(e,1),p.data), tper, expRate.xy);
%tind = cf(@(e,o) e(1:size(o,1),:), tind, obsCount);

ddzMask       = cf(@(e,o) e(1:size(o,1),:), ddzMask,       tind);
expCount.xy   = cf(@(e,o) e(1:size(o,1),:), expCount.xy,   tind);
expCount.xyhb = cf(@(e,o) e(1:size(o,1),:), expCount.xyhb, tind);

obsCount = cf(@(e) ...
              sq(sum(GetSegs(e.data,1:sampleRate.*ufrBinWidth:size(e,1),...
                             sampleRate.*ufrBinWidth,0),'omitnan')), ...
              ufr);


S = obsCount{t}(tind{t}&ddzMask{t}(:,u),u);
N = expCount.xy{t}(tind{t}&ddzMask{t}(:,u),u);        
N = expCount.xyhb{t}(tind{t}&ddzMask{t}(:,u),u);        
figure,plot(S,N,'.'),daspect([1,1,1])


%Z = (S-N+sum([-0.5,0.5].*double([S>=N,S<N])))/sqrt(N);
clear Z
Z.xy = [];
Z.xyhb = [];
ucnt = 0;
for t = [1:23];
    for u = 1:numel(units{t}),
        if pft{t}.maxRate(units{t}(u))>2;
            S = round(obsCount{t}(tind{t}&ddzMask{t}(:,u),u));
            N = expCount.xy{t}(tind{t}&ddzMask{t}(:,u),u);        
            Z.xy = [Z.xy;(S-N+sum(bsxfun(@times,[-0.5,0.5],double([S>=N,S<N])),2))./sqrt(N)];
            N = expCount.xyhb{t}(tind{t}&ddzMask{t}(:,u),u);        
            Z.xyhb = [Z.xyhb;(S-N+sum(bsxfun(@times,[-0.5,0.5],double([S>=N,S<N])),2))./sqrt(N)];
            ucnt = ucnt + 1;
        end
    end
end

figure();
subplot(211);
bar(linspace(-12,12,120),histc(Z.xy,  linspace(-12,12,120)),'histc');
hold('on');
plot(linspace(-12,12,120),1500*normpdf(linspace(-12,12,120),0,1),'m');
subplot(212);
bar(linspace(-12,12,120),histc(Z.xyhb,linspace(-12,12,120)),'histc');
hold('on');
plot(linspace(-12,12,120),1500*normpdf(linspace(-12,12,120),0,1),'m');