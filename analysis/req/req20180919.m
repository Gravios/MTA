% req20180919
%  Status: active
%  Type: Analysis
%  Author: Justin Graboski
%  Final_Forms: NA
%  Project: MjgER2016:placefields
%  Description: Theta phase precession of a place cells spiking activity seems to be independent of 
%               non spatial modulators of firing rate.
%  Protocol: 
%    1. select only first one or two spikes of each theta cycle.
%    2. compute placefields

MjgER2016_load_data;

states = {'theta-groom-sit','rear','hloc','hpause','lloc',...
          'lpause'};


sampleRate = [];

trialIndex = 20; % jg05-20120312.cof.all
Trial = Trials{trialIndex}; 
xyz = preproc_xyz(Trial,'trb');            
xyz.resample(sampleRate);

mazeCenterDist = sqrt(sum(xyz(:,'hcom',[1,2]).^2,3));
mazeCenterAng = circ_dist(atan2(xyz(:,'hcom',2),xyz(:,'hcom',1)),...
                atan2(diff(xyz(:,{'hcom','nose'},2),1,2),diff(xyz(:,{'hcom','nose'},1),1,2)));

stc = Trial.load('stc','msnn_ppsvd_raux');


unitSubset = units{trialIndex};
pft = pfs_2d_theta(Trial,unitSubset);
spk = Trial.load('spk',xyz.sampleRate,'gper',unitSubset,'deburst');
lfp = Trial.load('lfp',sessionList(trialIndex).thetaRef);
phz = phase(lfp,[5,12]);
phz.data = unwrap(phz.data);
phz.resample(xyz);
phzInd = discretize(phz.data,0:2*pi:max(phz.data(:)));
phz.data = mod(phz.data+pi,2*pi)-pi;
%lfp.resample(xyz);    



stcm = stc2mat(stc,xyz,states);

ddz            = compute_ddz(Trial,unitSubset,pft);%,[],[],'nose',[],copy(xyz),[],[],sampleRate);
[drz,~,drzAng] = compute_drz(Trial,unitSubset,pft);%,[],[],'nose',[],copy(xyz),[],[],sampleRate);
[hrz,~,hrzAng] = compute_hrz(Trial,unitSubset,pft);%,[],[],'nose',[],copy(xyz),[],[],sampleRate);

thresholds.mazeCenterDist = 400;
thresholds.mazeCenterAng = pi/2;
cind =    mazeCenterDist < thresholds.mazeCenterDist  |  abs(mazeCenterAng) < thresholds.mazeCenterAng ;


uind = find(unitSubset==21);%59,35
unit = unitSubset(uind);

res = spk(unit);
%res = res((stcm(res,1)==1)); % theta
res = res((stcm(res,3)==3|stcm(res,5)==5)&stcm(res,1)==1);
%res = res(stcm(res,3)==3&stcm(res,1)==1); % hloc
%res = res(stcm(res,5)==5&stcm(res,1)==1); % hloc
res(abs(ddz(res,uind))>300) = [];
if numel(res) > 10,
    res(res>xyz.size(1))=[];
    drzspk = drz(res,uind);
    hrzspk = hrz(res,uind);                
    ddzspk = ddz(res,uind);
    phzspk = phz(res,spk.map(spk.map(:,1)==unitSubset(uind),2));
    gind = ~isnan(drzspk)&~isnan(phzspk);                
    numRes = sum(gind);
else
    res = [];
    drzspk=[];
    ddzspk=[];
    phzspk=[];
    gind=[];
    numRes = 0;
end

[~,rind] = unique(phzInd(res,spk.map(spk.map(:,1)==unitSubset(uind),2)),'first');
%[~,rind] = unique(phzInd(res,spk.map(spk.map(:,1)==unitSubset(uind),2)),'last');
ures = res(rind);
if numel(res) >10,
    udrzspk = drz(ures,uind);
    uhrzspk = hrz(ures,uind);                
    uddzspk = ddz(ures,uind);
    uphzspk = phz(ures,spk.map(spk.map(:,1)==unitSubset(uind),2));
    ugind = ~isnan(udrzspk)&~isnan(uphzspk);
    numUres = sum(ugind);
else
    res = [];
    drzspk=[];
    ddzspk=[];
    phzspk=[];
    ugind=[];
    numUres = sum(ugind);    
end



figure,
subplot(221);
scatter([drzspk;drzspk],[phzspk;phzspk+2*pi],10,repmat(drzAng(res,uind),[2,1]),'filled');
subplot(222);
scatter([udrzspk;udrzspk],[uphzspk;uphzspk+2*pi],10,repmat(drzAng(ures,uind),[2,1]),'filled');
subplot(223);
scatter([hrzspk;hrzspk],[phzspk;phzspk+2*pi],10,repmat(hrzAng(res,uind),[2,1]),'filled');
subplot(224);
scatter([uhrzspk;uhrzspk],[uphzspk;uphzspk+2*pi],10,repmat(hrzAng(ures,uind),[2,1]),'filled');
colormap('hsv');

