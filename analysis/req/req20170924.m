% $$$ Trials = {};
% $$$ Trials{end+1} = 'jg05-20120309.cof.all';
% $$$ anatGrpSWPCenter{end+1}(1:5) = nan;
% $$$ anatGrpSWPCenter{end}(6) = 48;
% $$$ anatGrpSWPCenter{end}(7) = 55;
% $$$ anatGrpSWPCenter{end}(8) = 64;
% $$$ anatGrpSWPCenter{end}(9) = 70;
% $$$ anatGrpSWPCenter{end}(10) = 70;
% $$$ 
% $$$ Trials{end+1} = 'jg05-20120310.cof.all';
% $$$ anatGrpSWPCenter{end+1}(1:5) = nan;
% $$$ anatGrpSWPCenter{end}(6) = 48;
% $$$ anatGrpSWPCenter{end}(7) = 56;
% $$$ anatGrpSWPCenter{end}(8) = 63;
% $$$ anatGrpSWPCenter{end}(9) = 70;
% $$$ anatGrpSWPCenter{end}(10) = 70;
% $$$ 
% $$$ Trials{end+1} = 'jg05-20120311.cof.all';
% $$$ anatGrpSWPCenter{end+1}(1:3) = nan;
% $$$ anatGrpSWPCenter{end}(4) = 32;
% $$$ anatGrpSWPCenter{end}(5) = 39;
% $$$ anatGrpSWPCenter{end}(6) = 43;
% $$$ anatGrpSWPCenter{end}(7) = 49;
% $$$ anatGrpSWPCenter{end}(8) = 57;
% $$$ anatGrpSWPCenter{end}(9) = 69;
% $$$ anatGrpSWPCenter{end}(10) = 69;
% $$$ 
% $$$ Trials{end+1} = 'jg05-20120312.cof.all';
% $$$ anatGrpSWPCenter{end+1}(1:3) = nan;
% $$$ anatGrpSWPCenter{end}(4) = 32; % +2?
% $$$ anatGrpSWPCenter{end}(5) = 40; % +1?
% $$$ anatGrpSWPCenter{end}(6) = 46;
% $$$ anatGrpSWPCenter{end}(7) = 53;
% $$$ anatGrpSWPCenter{end}(8) = 61;
% $$$ anatGrpSWPCenter{end}(9) = 71;
% $$$ anatGrpSWPCenter{end}(10) = 71;
% $$$ 
% $$$ Trials{end+1} = 'jg05-20120317.cof.all';
% $$$ anatGrpSWPCenter{end+1}(1:5) = nan;
% $$$ anatGrpSWPCenter{end}(6) = 46;
% $$$ anatGrpSWPCenter{end}(7) = 51;
% $$$ anatGrpSWPCenter{end}(8) = 57;
% $$$ anatGrpSWPCenter{end}(9) = 70;
% $$$ anatGrpSWPCenter{end}(10) = 70;

Trials = {};
anatGrpSWPCenter = {};
Trials{end+1} = 'jg05-20120309.cof.all';
anatGrpSWPCenter{end+1}(1:5) = nan;
anatGrpSWPCenter{end}(6) = 8;
anatGrpSWPCenter{end}(7) = 7;
anatGrpSWPCenter{end}(8) = 8;
anatGrpSWPCenter{end}(9) = 7;
anatGrpSWPCenter{end}(10) = 1;
anatGrpSWPCenter{end}(11:13) = nan;

Trials{end+1} = 'jg05-20120310.cof.all';
anatGrpSWPCenter{end+1}(1:5) = nan;
anatGrpSWPCenter{end}(6) = 8;
anatGrpSWPCenter{end}(7) = 8;
anatGrpSWPCenter{end}(8) = 7;
anatGrpSWPCenter{end}(9) = 7;
anatGrpSWPCenter{end}(10) = 3;
anatGrpSWPCenter{end}(11:12) = nan;

Trials{end+1} = 'jg05-20120311.cof.all';
anatGrpSWPCenter{end+1}(1:3) = nan;
anatGrpSWPCenter{end}(4) = 8;
anatGrpSWPCenter{end}(5) = 7;
anatGrpSWPCenter{end}(6) = 3;
anatGrpSWPCenter{end}(7) = 1;
anatGrpSWPCenter{end}(8) = 1;
anatGrpSWPCenter{end}(9) = 6;
anatGrpSWPCenter{end}(10) = 2;
anatGrpSWPCenter{end}(11:12) = nan;

Trials{end+1} = 'jg05-20120312.cof.all';
anatGrpSWPCenter{end+1}(1:3) = nan;
anatGrpSWPCenter{end}(4) = 10; % +2 greater than chan count estimate 
anatGrpSWPCenter{end}(5) = 9;  % +1 greater than chan count estimate  
anatGrpSWPCenter{end}(6) = 8;
anatGrpSWPCenter{end}(7) = 5;
anatGrpSWPCenter{end}(8) = 5;
anatGrpSWPCenter{end}(9) = 8;
anatGrpSWPCenter{end}(10) = 4;
anatGrpSWPCenter{end}(11:12) = nan;

Trials{end+1} = 'jg05-20120317.cof.all';
anatGrpSWPCenter{end+1}(1:5) = nan;
anatGrpSWPCenter{end}(6) = 6;
anatGrpSWPCenter{end}(7) = 3;
anatGrpSWPCenter{end}(8) = 1;
anatGrpSWPCenter{end}(9) = 7;
anatGrpSWPCenter{end}(10) = 3;
anatGrpSWPCenter{end}(11:12) = nan;

Trials = cf(@(t)    MTATrial.validate(t),       Trials);
cf(@(t)    t.load('stc','msnn_ppsvd'), Trials);
%cf(@(t)    t.load('stc','NN0317R'), Trials);
pft    = cf(@(t)    pfs_2d_theta(t,[],[],false), Trials);
mrt    = cf(@(p)    p.maxRate,                  pft);
units =  cf(@(t)    select_units(t,18),         Trials);
units =  cf(@(u,m,p)  u(m(p.data.clu(u))>2),    units, mrt, pft);

states = {'walk&theta','rear&theta'};
%states = {'lloc&theta','hloc&theta','rear&theta'};

% Compute place fields and subsampled estimate
pfkbs = {};
for t = 1:numel(Trials),
    for sts = 1:numel(states),    
        defargs = get_default_args_MjgEdER2016('MTAApfs','struct');
        defargs.units = units{t};
        defargs.states = states{sts};
        %defargs.overwrite = true;
        defargs = struct2varargin(defargs);
        pfkbs{t,sts} = MTAApfs(Trials{t},defargs{:});      
    end
end



% UPDATED MTA/utilities/NeuronQuality.m
% added field to output which indicates the channel where a unit's peak amplitude occured

%nq = cf(@(t) NeuronQuality(t,[],[],[],false),  Trials);
cf(@(t) t.load('nq'), Trials);
nq = cf(@(t) t.nq, Trials);

allUnitDepths = cf(@(t,n,a)  n.maxAmpChan-a(n.ElNum)',  Trials,nq,anatGrpSWPCenter);


maxUnitRates = cf(@(p)    p.maxRate,                  pfkbs);
meanUnitRates = cf(@(w,r)  (w+r)./2,         maxUnitRates(:,1),maxUnitRates(:,2));
rateRatio    = cf(@(w,r)  (w-r)./(w+r),      maxUnitRates(:,1),maxUnitRates(:,2));
unitDepths   = cf(@(a,u)  a(u),                       allUnitDepths,units);
elnum        = cf(@(p)    p.data.el',                  pfkbs(:,1));
%spInfo       = cf(@(p)    p.data.el',                 pft;

en = cat(1,elnum{:});
mur = cat(1,meanUnitRates{:});
rr = cat(1,rateRatio{:});
ud = cat(1,unitDepths{:});

ud(en<9) = ud(en<9)*20;
ud(en>9) = ud(en>9)*50;

ind = mur>8;
figure,hold on
plot(rr(ind)+randn(size(rr(ind)))/20,ud(ind)+randn(size(ud(ind)))*2,'b.')


figure();
hold('on');
plot(log10(mur(ind))+randn(size(mur(ind)))/10,ud(ind)+randn(size(ud(ind)))*2,'.b')


figure();
hold('on');
plot(log10(mur(ind))+randn(size(mur(ind)))/10,rr(ind)+randn(size(rr(ind)))/20,'.b')
