MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
 
Trial = MTATrial('jg05-20120317','all');
Trial.xyz.load(Trial);
Trial.filter('xyz',gausswin(31)./sum(gausswin(31)));





%% Figure 1 3D tracking -> labeled behavior

figure,
sph = [];
% A: Image, rat with markers in maze.
%subplot2(6,6,1,[1,2]); image('rat.png');
% B: Image, marker skeleton reconstruction
%subplot2(6,6,1,[3,4]); image('skel.png');
% C: Image, cameras and maze
%subplot2(6,6,1,[3,4]); image('skel.png');

% D: Plot ( Time Series ), COM_head speed ( blue )
%    Plot ( Time Series ), COM_body speed ( green )
%    Plot ( Time Series ), COM_head height ( red )

window_figD = 1300; index_figD = 600;

coms = cat(2,Trial.com(Trial.xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})),...
             Trial.com(Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'})));
t = linspace(0,size(coms,1)/Trial.xyz.sampleRate,size(coms,1)-1)';
vcoms = sqrt(sum(diff(coms).^2,3));
sph(end+1) = subplot2(6,6,2,[1:6]); 
ax_figD = plotyy(t,clip(coms(1:end-1,2,3),20,350),[t,t],clip(vcoms,0,25));
set(ax_figD,'xlim',[479,515]);%[1660,1730]);
set(ax_figD(1),'ylim',[20,280]);
set(ax_figD(2),'ylim',[-2,10]);



%% Figure 1 BHV state Place Fields and ccgs
% rear walk hwalk lwalk

Trial = MTATrial('jg05-20120309','all');
Trial.load('nq');
units = find(Trial.nq.SpkWidthR>0.8&Trial.nq.eDist>18)';
states = {'rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
nsts = numel(states);
for i = 1:nsts,
pfs{i} =     MTAAknnpfs(Trial,units,states{i},0,'numIter',1000, ...
                'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
end


f = figure(1123);
for u = pfs{1}.data.clu, 
for i = 1:nsts,
subplotfit(i,nsts),pfs{i}.plot(u);title([pfs{i}.parameters.states,':',num2str(u)]),
end
pause(.2),
reportfig(f,[Trial.name,'-pfs_state'],[],['unit: ',num2str(u)]);
end


figc = figure(1231);
for u = pfl.data.clu, 
for i = 1:nsts,
subplotfit(i,nsts),pfs{i}.plot(u);title([pfs{i}.parameters.states,':',num2str(u)]),
end
pause(.2),
reportfig(figc,'jg05-201203017-ccg_state',[],['unit: num2str(u)']);
end

bccg{i} = gen_bhv_ccg(Trial);
figure,
subplot(211),bccg.plot(units(10),1);
subplot(212),bccg.plot(units(10),2);


%% End - Figure 1




%% Figure 2 - Marker vs Marker Speed Joint Distribution
MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317','all');
Trial.xyz.load(Trial);
Trial.filter('xyz',gausswin(31)./sum(gausswin(31)));

states = 'twr';
states = 'wgl';
m1 = 'spine_lower';
m2 = 'head_front';


nbins=100;mab = 2.1;mib = -2;
nbins=50;mab = 2.1;mib = 0.5;

hc = cat(1,mib,mab,mib,mab);
bc = cat(1,mab,mib,mib,mab);


dmax = 100;dmin =0.1;
dxys = Trial.vel(1:Trial.xyz.size(2),[1,2]).*Trial.xyz.sampleRate./10;
dxys = MTADxyz([],[],dxys,Trial.xyz.sampleRate,[],[],Trial.xyz.model);

figure,set(0,'defaulttextinterpreter','none')
for state=states,subplot(3,1,find(state==states));

b = clip(log10(dxys(Trial.stc{state},m1)),mib,mab);
h = clip(log10(dxys(Trial.stc{state},m2)),mib,mab);

hist2([[b;bc],[h;hc]],nbins,nbins);
title([m1 ' vs ' m2 ' during ' Trial.stc{state}.label])
xlabel([m1 ' speed (cm/s)'])
ylabel([m2 ' speed (cm/s)'])
caxis([0,250]);
line([mib,mab],[mib,mab],'color',[1,1,1]);

xticks = get(gca,'XTickLabel');
xticks = mat2cell(xticks,ones(1,size(xticks,1)),size(xticks,2));
xticks = cellfun(@str2num,xticks);
xticks = mat2cell(10.^xticks,ones(1,size(xticks,1)),1);
xticks = cellfun(@sprintf,repmat({'%2.2f'},numel(xticks),1),xticks,'uniformoutput',false);
set(gca,'XTickLabel',xticks);

yticks = get(gca,'YTickLabel');
yticks = mat2cell(yticks,ones(1,size(yticks,1)),size(yticks,2));
yticks = cellfun(@str2num,yticks);
yticks = mat2cell(10.^yticks,ones(1,size(yticks,1)),1);
yticks = cellfun(@sprintf,repmat({'%2.2f'},numel(yticks),1),yticks,'uniformoutput',false);
set(gca,'YTickLabel',yticks);

end


%% End - Figure 2

%% Figure 3 - PlaceFields of behaviors for multiple units sorted by behavioral modulation

MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
sname = 'jg05-20120310';
Trial = MTATrial(sname,'all');
Trial.load('nq');

states = {'walk&theta','hwalk&theta','lwalk&theta'};
%states = {'walk&theta','hwalk&theta','lwalk&theta'};

% Get the place fields
pfs = {};
for s = 1:numel(states);
    pfs{s} = MTAAknnpfs(Trial,[],states{s},0,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThresh',70);
end

% Get the auto correlograms
[accg,tbins] = autoccg(MTASession(sname));

% Select Some units
spc = [];
for i=1:numel(pfs),
    spc = cat(2,spc,pfs{i}.spatialCoherence([]));
end
mxr = [];
for i=1:numel(pfs),
    mxr = cat(2,mxr,pfs{i}.maxRate([]));
end

ind = (spc(:,2)>.98|spc(:,3)>.98)&Trial.nq.eDist>30&Trial.nq.SpkWidthR>.5;
figure,plot(log10(mxr(ind,2)),log10(mxr(ind,3)),'.')
line([-1,2],[-1,2]);

wunits = find(ind&log10(mxr(:,2))<0);
wrunits = find(ind&log10(mxr(:,2))>.8&log10(mxr(:,3))>.8);

[~,wrind] =sort(mxr(wrunits,2)./mxr(wrunits,3));
wrunits(wrind) = wrunits;

runits = find(ind&log10(mxr(:,3))<.54);

sunits = [wunits;wrunits;runits];

% Plot it all
nrows = numel(states)+1;
ncols = numel(sunits);

figure
for unit = sunits'
    un = find(unit==sunits);
    mufr = max(mxr(unit,:))
    for s = 1:nrows-1;
        subplot2(nrows,ncols,s,un);
        pfs{s}.plot(unit,[],[],[0,mufr]);
    end
    subplot2(nrows,ncols,s+1,un);
    bar(tbins,accg(:,unit));axis tight;
end


%% End Figure - 3



%% Figure 4 - Classical 2-D phase precession
%sname = 'jg05-20120317';
%sname = 'jg05-20120310';
sname = 'er06-20130612';
%sname = 'jg05-20120309';

Trial = MTATrial(sname,'all');
Trial.xyz.load(Trial);

%Trial.load('nq');
%Trial.lfp.load(Trial,[61,71:3:84]);
Trial.lfp.load(Trial,[32]);



tbp_phase = Trial.lfp.phase;


state = 'walk';

tbp_phase.resample(Trial.xyz);


pfs = MTAAknnpfs(Trial,[],state,0,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThresh',100,'nNearestNeighbors',120);

Trial.spk.create(Trial,Trial.xyz.sampleRate,state,[],'deburst');

units = pfs.data.clu(:)';

% Get the expected ufr for each xy 
wpmr = zeros(Trial.xyz.size(1),numel(units));
wpmr = ones(Trial.xyz.size(1),numel(units));
[~,indx] = min(abs(repmat(pfs.adata.bins{1}',Trial.xyz.size(1),1)-repmat(Trial.xyz(:,7,1),1,numel(pfs.adata.bins{1}))),[],2);
[~,indy] = min(abs(repmat(pfs.adata.bins{2}',Trial.xyz.size(1),1)-repmat(Trial.xyz(:,7,2),1,numel(pfs.adata.bins{2}))),[],2);
indrm = sub2ind(pfs.adata.binSizes',indx,indy);


for unit = units,
    %rateMap = reshape(pfs.data.rateMap(:,unit),pfs.adata.binSizes');
    rateMap = pfs.plot(unit);
    %rateMap = rot90(rot90(rateMap)',1);
    wpmr(:,unit==units) = rateMap(indrm);
end

%wpmr(wpmr<1)=1;
%figure,scatter(Trial.xyz(1:61:end,7,1),Trial.xyz(1:61:end,7,2),wpmr(1:61:end,20))

[pmr,pmp] = pfs.maxRate([]);
pmr = repmat(pmr(:)',Trial.xyz.size(1),1);
pmp = fliplr(pmp);

for unit = units
    %fprintf(['head pfs angle for unit: %i \n'],unit)
    %pfang = MTADang;
    pfhxy = Trial.xyz(:,{'head_back','head_front'},:);
    pfhxy = cat(2,pfhxy,permute(repmat([pmp(unit==units,:),0],Trial.xyz.size(1),1),[1,3,2]));
    pfhxy = MTADxyz([],[],pfhxy,Trial.xyz.sampleRate);

cor = cell(1,3);
[cor{:}] = cart2sph(pfhxy(:,2,1)-pfhxy(:,1,1),pfhxy(:,2,2)-pfhxy(:,1,2),pfhxy(:,2,3)-pfhxy(:,1,3));
cor = cell2mat(cor);

por = cell(1,3);
[por{:}] = cart2sph(pfhxy(:,3,1)-pfhxy(:,1,1),pfhxy(:,3,2)-pfhxy(:,1,2),pfhxy(:,3,3)-pfhxy(:,1,3));
por = cell2mat(por);

pfds(:,unit==units) = circ_dist(cor(:,1),por(:,1));

    %pfang.create(Trial,pfhxy);
    %pfd(:,unit==units) = circ_dist(pfang(:,1,2,1),pfang(:,1,3,1));
end

pfd = zeros(size(pfds));
pfd(abs(pfds)<=pi/2)=-1;
pfd(abs(pfds)>pi/2)=1;


%DRZ 
DRZ = pfd.*(1-wpmr./pmr);
%DRZ = pfd.*(1-log10(clip(wpmr,1,100))./log10(pmr));
%figure,scatter(Trial.xyz(1:31:end,7,1),Trial.xyz(1:31:end,7,2),clip(1./abs(DRZ(1:31:end,2)),0,50))

phchan = 1;

figure,
for unit = units
clf
res = Trial.spk(unit);
rind = SplitIntoBursts(res,3);
res = res(unique(rind));

if isempty(res),continue,end
if numel(res)>20000||numel(res)<100,continue,end

drzspk = DRZ(res,unit);
phzspk = tbp_phase(res,phchan);

gind = ~isnan(drzspk)&~isnan(phzspk);

subplot2(2,3,[1,2],1),
plot(Trial.xyz(res,7,1),Trial.xyz(res,7,2),'.');
xlim([-500,500]),ylim([-500,500])

subplot2(2,3,[1,2],2),pfs.plot(unit);
hold on,plot(pmp(unit==units,1),pmp(unit==units,2),'w*')
title(num2str(unit))

subplot2(2,3,1,3),plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'.')
hold on,          plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'.')
xlim([-1,1]),

subplot2(2,3,2,3)
hist2([[drzspk(gind);drzspk(gind)],[circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360]],30,25)



pause(.1)
waitforbuttonpress

reportfig(gcf,'er06-20130613-2Dphspredb',0,num2str(unit),[],0);

end


%% End Figure 4





%% Figure 5 - State Wise Unit Firing Rates
sname = 'jg05-20120309';
Trial = MTATrial(sname,'all');
Trial.load('nq');

states = {'theta','vel&theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
units = find(Trial.nq.eDist>30&Trial.nq.SpkWidthR>.5)';
sscount = nan(numel(units),numel(states));
ssdur   = nan(numel(units),numel(states));

for s = 1:numel(states),

    tsts = Trial.stc{states{s}};
    Trial.spk.create(Trial,Trial.xyz.sampleRate,states{s},units);
    ssdur(:,s) = sum(diff(tsts.data,1,2));
for u = units,

    try, sscount(u==units,s) = numel(Trial.spk(u));end

end
end
srates = sscount./(ssdur./Trial.xyz.sampleRate);

nrates = srates./repmat(max(srates,[],2),1,numel(states));

[~,rind] = sort(nrates(:,3));

figure,imagescnan(nrates(rind,:)',[],[],1);
set(gca,'YTickLabel',states);

%% End Figure 5




%% Figure 6 - 

MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%MTAConfiguration('/data/data/gravio','absolute');
sname = 'jg05-20120317';
%sname = 'jg05-20120310';
%sname = 'jg05-20120309';
disp = 0;
%chans = [1:2:8];
chans = [71:2:96];
marker = 'spine_lower'

Trial = MTATrial(sname,'all');
Trial.ang.load(Trial);
Trial.xyz.load(Trial);
Trial.lfp.load(Trial,chans);


wlfp = WhitenSignal(Trial.lfp.data,[],1);


matlabpool open 12

tl=[];
fl=[];
yl=[];
parfor i = 1:Trial.lfp.size(2),
    [yl(:,:,i),fl(:,i),tl(:,i)] = mtchglong(wlfp(:,i),2^12,Trial.lfp.sampleRate,2^11,2^11*0.875,[],[],[],[1,40]);
end

th=[];
fh=[];
yh=[];
parfor i = 1:Trial.lfp.size(2),
    [yh(:,:,i),fh(:,i),th(:,i)] = mtchglong(wlfp(:,i),2^9,Trial.lfp.sampleRate,2^8,2^8*0.875,[],[],[],[40,120]);
end
yld = MTADlfp('data',yl,'sampleRate',1/diff(tl(1:2,1)));
yhd = MTADlfp('data',yh,'sampleRate',1/diff(th(1:2,1)));

bang = ButFilter(Trial.ang(:,4,5,3),3,[1,20]./(Trial.ang.sampleRate./2),'bandpass');
bhh = ButFilter(Trial.xyz(:,7,3),3,[1,20]./(Trial.xyz.sampleRate./2),'bandpass');
bhx = ButFilter(Trial.xyz(:,7,1),3,[1,20]./(Trial.xyz.sampleRate./2),'bandpass');
bhy = ButFilter(Trial.xyz(:,7,2),3,[1,20]./(Trial.xyz.sampleRate./2),'bandpass');
if disp, figure,plot([bang,bhh,bhx,bhy]+3.*repmat(1:4,size(bang,1),1)), end

wang = WhitenSignal([bang,bhh,bhx,bhy],[],1);
[ya,fa,ta,phia,fsta] = mtchglong(wang,2^9,Trial.ang.sampleRate,2^8,2^8*0.875,[],[],[],[2,16]);

%figure,plot(mean(ya(:,fa>9&fa<12),2)./mean(ya(:,fa<7),2))
%figure,plot([mean(ya(:,fa>9&fa<12,1,1),2),mean(ya(:,fa>9&fa<12,2,2),2)])
%figure,plot(mean(ya(:,fa>9&fa<12,1,1),2),mean([mean(ya(:,fa>9&fa<12,4,4),2),mean(ya(:,fa>9&fa<12,3,3),2)],2),'.')
%figure,plot(log10(mean(ya(:,fa>9&fa<12,1,1),2)),log10(mean([mean(ya(:,fa>9&fa<12,4,4),2),mean(ya(:,fa>9&fa<12,3,3),2)],2)),'.')
spowa = log10(mean(ya(:,fa>6&fa<12,1,1),2));
spowa = log10(mean(ya(:,fa>6&fa<12,2,2),2));



% VELOCITY 
xyz = Trial.xyz.copy;
xyz.filter(gausswin(31)./sum(gausswin(31)));
v = MTADxyz([],[],sqrt(sum(diff(xyz(:,marker,[1,2])).^2,3)).*Trial.xyz.sampleRate./10,Trial.xyz.sampleRate);



%figure,hist2([clip(log10(v.data),-2,3),clip(spow,-6,1)],100,100),caxis([0,40])
%figure,hist2([clip(log10(v.data),-2,3),clip(mean([spowa,spowh],2),-6,1)],100,100),caxis([0,40])

yad = MTADlfp([],[],ya,1/diff(ta(1:2)));

% RESAMPLE variables
xyz.resample(yad);
v.resample(yad);
yld.resample(yad);
%hdv.resample(yad);


sbins = 25;
sedges = [-5,-2];
%sedges = [50,160];

vbins = 25;
vedges = [0,2];

wper = Trial.stc{'w'}.copy;
wper.cast('TimeSeries');
wper.resample(yad);

%spow = xyz(:,7,3);
spow = spowa;
spow = clip(spow,sedges(1),sedges(2));
sind = spow>sedges(1)&spow<sedges(2);
vlog = clip(log10(v.data),vedges(1),vedges(2));
vind = vlog>vedges(1)&vlog<vedges(2)&~isnan(vlog);
aind = sind&vind&wper.data;

[spow_count,shind] = histc(spow(aind),sedges(1):abs(diff(sedges))/sbins:sedges(2));
[v_count,vhind] = histc(vlog(aind),vedges(1):abs(diff(vedges))/vbins:vedges(2));

%tpow = log10(mean(yld(aind,fl>=6&fl<=12,3),2));
tpow = log10(mean(yld(aind,fl>=4&fl<=16,1,1),2));
%tpow = log10(mean(yad(aind,fa>=4&fa<=16,1,1),2));
tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;

%A = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@mean,nan);
%figure,
%subplot(133),imagescnan({vedges,sedges,A'},[-5,-3],[],1,[0,0,0]),axis xy,
%clf
%imagescnan({vedges,sedges,A'},[1.2,1.9],[],1,[0,0,0]),axis xy,
%subplot(122),imagescnan({vedges,sedges,A'},[],[],1,[0,0,0]),axis xy,

% CA1 LM 81;
% DG  G  85;
% CA3 ?  95;
chan = find(chans == 71);
numIter = 10000;
%tpow = log10(mean(yld(aind,fl>6&fl<12,chan),2));
%tpow = log10(mean(yld(aind,fl<=4,chan),2));
tpow = log10(mean(yad(aind,fa>=4&fa<=16,1,1),2));
%tpow = log10(mean(yld(aind,fh>50&fh<80,chan),2));

B=nan(vbins,sbins,numIter);
A=nan(vbins,sbins,numIter);
%S=nan(vbins,sbins,numIter);
tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
B = accumarray([vhind(tind),shind(tind)],ones(sum(tind),1),[vbins,sbins],@sum,nan);
A(:,:,1) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@median,nan);
%S(:,:,1) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@std,nan);
for i = 2:numIter,
tpow = tpow(randperm(numel(tpow)));
tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
A(:,:,i) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@median,nan);
%S(:,:,i) = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@std,nan);
end

AS = sort(A,3);
P = 1./sum(repmat(A(:,:,1),[1,1,numIter])>A,3);
P(isinf(P)) = nan;

SIG = P<=0.0002;
ASIG = A; 
ASIG(~SIG)=nan;
ASIG(B<10)=nan;
Aclims = [prctile(ASIG(~isnan(ASIG)),5),prctile(ASIG(~isnan(ASIG)),95)];

figure

subplot(131),imagescnan({vedges,sppedges,B'./yad.sampleRate},[],[],1,[0,0,0]),axis xy,
title('Occupancy in seconds')
ylabel('log10(P(10Hz)) of Spine to Head Distance')
xlabel('Head Speed (cm/s)')
ticks_lin2log(gca,'x')

subplot(132),imagescnan({vedges,sedges,   A(:,:,1)'},Aclims,[],1,[0,0,0]),axis xy,
title('Mean Power 1-4Hz given 10Hz Osc. Power dist(SU,HB) VS Vel(HF)')
ylabel('log10(P(10Hz)) of Spine to Head Distance')
xlabel('Head Speed (cm/s)')
ticks_lin2log(gca,'x')

subplot(133),imagescnan({vedges,sedges,ASIG(:,:,1)'},Aclims,[],1,[0,0,0]),axis xy,
title('P<0.0002 and Occupancy > 1.33 seconds')
ylabel('log10(P(10Hz)) of Spine to Head Distance')
xlabel('Head Speed (cm/s)')
ticks_lin2log(gca,'x')





flim=[1,4;6,12;20,27;30,40];
%flim=[40,60;60,80;80,100;100,120];
%mychans = [71,73,81,85,95];
mychans = 1:2:8;
B = accumarray([vhind(tind),shind(tind)],ones(sum(tind),1),[vbins,sbins],@sum,nan);
figure
for c = 1:numel(mychans),
    for i = flim',
        subplot2(numel(mychans),size(flim,1),c,find(i(1)==flim(:,1)));
        tpow = log10(mean(yld(aind,fl>i(1)&fl<i(2),find(chans==mychans(c))),2));
        %tpow = log10(mean(yld(aind,fh>i(1)&fh<i(2),find(chans==mychans(c))),2));
        tind = ~isinf(tpow)&~isnan(tpow)&vhind~=0&shind~=0;
        AFB = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@mean,nan);
        AFB(B<5)=nan;
        AFBclims = [prctile(AFB(~isnan(AFB)),5),prctile(AFB(~isnan(AFB)),95)];
        imagescnan({vedges,sedges,AFB'},AFBclims,[],1,[0,0,0]);
        axis xy,
        ticks_lin2log(gca,'x')
        title(['C: ' num2str(mychans(c)) 'Mean P(' num2str(i(1)) '-' num2str(i(2)) ')'])
    end
end

text(.1,.1,['Mean P(' num2str(flim(1)) '-' num2str(flim(end)) ') given 10Hz Osc. Power dist(SU,HB) VS Vel(SL)'])
ylabel('log10(P(10Hz)) of Spine to Head Distance')
xlabel('Head Speed (cm/s)')


A = accumarray([vhind(tind),shind(tind)],tpow(tind),[vbins,sbins],@nanmean,nan);
figure,imagescnan({vedges,sedges,clip(A,0,.015)'},[],[],1,[0,0,0]),axis xy,
% $$$ 
% $$$ %tbp_phase.resample(yad);
% $$$ A = accumarray([vhind(tind),shind(tind)],tbp_phase(tind),[vbins,sbins],@circ_median,nan);
% $$$ A = accumarray([vhind(tind),shind(tind)],tbp_phase(tind),[vbins,sbins],@circ_var,nan);
% $$$ figure,imagescnan({vedges,sedges,A'},[],[],1,[0,0,0]),axis xy,


% $$$ %dratio = mean(yld(aind,fl>6&fl<12,chan),2)./mean(yld(aind,fl<4|(fl>12&fl<18),chan),2);
% $$$ %dratio = mean(yld(aind,fl>6&fl<12,chan),2)./mean(yld(aind,fl<4,chan),2);



%% End - Figure 6


%% Figure 7 - rhm vs thetaM
MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%MTAConfiguration('/data/data/gravio','absolute');
%sname = 'jg05-20120317';
sname = 'jg05-20120310';
%sname = 'jg05-20120309';
chans = [71:2:96];
marker = 'spine_lower'

Trial = MTATrial(sname,'all');
Trial.ang.load(Trial);
Trial.xyz.load(Trial);
Trial.lfp.load(Trial,chans);

wlfp = WhitenSignal(Trial.lfp.data,[],1);

yl=[];
for i = 1:Trial.lfp.size(2),
    [yl(:,:,i),fl,tl] = mtchglong(wlfp(:,i),2^12,Trial.lfp.sampleRate,2^11,2^11*0.875,[],[],[],[1,40]);
end
%yh=[];
%for i = 1:Trial.lfp.size(2),
%    [yh(:,:,i),fh,th] = mtchglong(wlfp(:,i),2^9,Trial.lfp.sampleRate,2^8,2^8*0.875,[],[],[],[40,120]);
%end

bang = ButFilter(Trial.ang(:,4,5,3),3,[1,20]./(Trial.ang.sampleRate./2),'bandpass');
%bhh = ButFilter(Trial.xyz(:,7,3),3,[2,16]./(Trial.xyz.sampleRate./2),'bandpass');
%bhx = ButFilter(Trial.xyz(:,7,1),3,[2,16]./(Trial.xyz.sampleRate./2),'bandpass');
%bhy = ButFilter(Trial.xyz(:,7,2),3,[2,16]./(Trial.xyz.sampleRate./2),'bandpass');
%if disp, figure,plot([bang,bhh,bhx,bhy]+3.*repmat(1:4,size(bang,1),1)), end

%wang = WhitenSignal([bang,bhh,bhx,bhy],[],1);
wang = WhitenSignal([bang],[],1);
[ya,fa,ta,phia,fsta] = mtchglong(wang,2^9,Trial.ang.sampleRate,2^8,2^8*0.875,[],[],[],[2,16]);
yad = MTADlfp([],[],ya,1/diff(ta(1:2)));
yld = MTADlfp([],[],yl,1/diff(tl(1:2)));
yld.resample(yad);
xyz = Trial.xyz.copy;
xyz.filter(gausswin(31)./sum(gausswin(31)));
v = MTADxyz([],[],sqrt(sum(diff(xyz(:,marker,[1,2])).^2,3)).*Trial.xyz.sampleRate./10,Trial.xyz.sampleRate);
v.resample(yad);

count = 1;
states = 'tvrwgl';
for c = 1:13,
for s = 1:numel(states)

wper = Trial.stc{states(s)}.copy;
wper.cast('TimeSeries');
wper.resample(yad);



ylpf = yld(:,fl<14&fl>5,c);yapf = yad(:,fa<14&fa>5,1,1);
%yad(:,fa<14&fa>5,3,3)+yad(:,fa<14&fa>5,4,4)+yad(:,fa<14&fa>5,1,1);
flpf = fl(fl<14&fl>5);fapf = fa(fa<14&fa>5);
[ylpfsp,ylpfs] = sort(ylpf,2,'descend');[yapfsp,yapfs] = sort(yapf,2,'descend');
flf = flpf(ylpfs(:,1));faf = fapf(yapfs(:,1));

aind = wper.data==1;
%aind = true(size(wper.data));
lbounds = prctile(log10(ylpfsp(:,1)),[.5,99.5]);
abounds = prctile(log10(yapfsp(:,1)),[.5,99.5]);
aind = aind&~isinf(log10(ylpfsp(:,1)))&~isinf(log10(yapfsp(:,1)));
%aind = aind&~isinf(log10(ylpfsp(:,1)))&~isinf(log10(yapfsp(:,1)))...
%        &log10(ylpfsp(:,1))>lbounds(1)&log10(ylpfsp(:,1))<lbounds(2)...
%        &log10(yapfsp(:,1))>abounds(1)&log10(yapfsp(:,1))<abounds(2);

%figure,hist(flf,29)
%figure,hist(faf,38)
%figure,hist2([flf(aind),faf(aind)],29,38);caxis([0,140])

%figure,hist2([clip(log10(ylpfsp(aind,1)),2.4,4), ...
%              clip(log10(yapfsp(aind,1)),-5.5,-2.5)],30,30);caxis([0,130])

%subplot(13,6,count),
%figure,
%plot(log10(ylpfsp(aind,1)),log10(yapfsp(aind,1)),'.')
for shift = 1:50;
[rho(c,s,shift),p(c,shift)] = corr(circshift(log10(ylpfsp(aind,1)),25-shift),log10(yapfsp(aind,1)));
end
%legend(num2str([rho,p]))
%count = count+1;
end
end
figure,imagesc(rho)


figure,
for i= 1:25
subplot(5,5,i),hist2([circshift(flf,-i),faf],29,38);caxis([0,150])
end

%% End - Figure 7

%% Figure 8 Feature Selection 
MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%MTAConfiguration('/data/data/gravio','absolute');
%sname = 'jg05-20120317';
sname = 'jg05-20120310';
%sname = 'jg05-20120309';
disp = 0;
Trial = MTATrial(sname);
Trial.load('xyz');
Trial.load('ang');

Trial.xyz.filter(gausswin(61)./sum(gausswin(61)));
vel = Trial.vel;
vel = clip(log10(vel),-2,2);
vel =  MTADxyz('data',[zeros([1,size(vel,2)]);vel],'sampleRate',Trial.xyz.sampleRate);

m=1;
figure
%hist2([vel(:,m),vel(:,7)],linspace(-2,2,64),linspace(-2,2,64));
hist2([vel(vel(:,1)>0,m),vel(vel(:,1)>0,7)],linspace(-2,2,64),linspace(-2,2,64));
caxis([0,600])


% marker speed vs marker height
figure
ind = ':';
vm = 1;
hm = 1;
hist2([vel(:,vm),clip(log10(Trial.xyz(:,hm,3)),0,3)],linspace(-2,2,64),linspace(1,2,64))
caxis([0,1600])
hold on
states = 'rwgl';
colors = 'wgmc';
for i = 1:numel(states),
ind = Trial.stc{states(i)};
vh = [vel(ind,vm),clip(log10(Trial.xyz(ind,hm,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style',colors(i));
end


% marker speed vs marker segment pitch
figure
ind = ':';
vm = 1;
ms = [3,4];
hist2([vel(:,vm),Trial.ang(:,ms(1),ms(2),2)],linspace(-1.5,2,64),linspace(-1.8,1.8,64))
%ticks_lin2log
caxis([0,1600])
hold on
states = 'rwgl';
colors = 'wgmc';
for i = 1:numel(states),
ind = Trial.stc{states(i)};
vh = [vel(ind,vm),Trial.ang(ind,ms(1),ms(2),2)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style',colors(i));
end



%% End - Figure 8

%% Figure - 9 Comodugram rhm and lfp


MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%MTAConfiguration('/data/data/gravio','absolute');

%sname = 'jg05-20120315';
sname = 'jg05-20120310';
%sname = 'jg05-20120309';
%sname = 'jg04-20120129';
%sname = 'jg04-20120130';
%sname = 'co01-20140222';

%sname = 'jg05-20120317';
chans = 68:3:95;

Trial = MTATrial(sname,'all');
%Trial.ang.load(Trial);
Trial.xyz.load(Trial);
Trial.lfp.load(Trial,chans);
Trial.lfp.resample(Trial.xyz);

rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = Trial.com(rb);
Trial.addMarker(Trial.xyz,'hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
Trial.addMarker(Trial.xyz,'fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(gausswin(61)./sum(gausswin(61)),hcom),[1,3,2]));

ang = Trial.ang.copy;
ang.create(Trial);
ang.data = ang(:,7,11,3);
%ang.resample(Trial.lfp);
ang = ang.data;
ang(isnan(ang))=nanmean(ang);

%bang = ButFilter(ang,3,[1,20]./(Trial.ang.sampleRate./2),'bandpass').*500;
bang = ang;

%abang = WhitenSignal(bang);
%lbang = WhitenSignal(Trial.lfp.data,[],1);
lbang = WhitenSignal([Trial.lfp.data,bang],[],1);
%lbang = MTADlfp([],[],lbang,Trial.lfp.sampleRate);
lbang = MTADlfp([],[],lbang,Trial.ang.sampleRate);
%lbang = MTADlfp([],[],[lbang,abang],Trial.ang.sampleRate);



%[ta,fa,ya,pha,sfa] = mtchglong(lbang(:,[1,9]),2^9,Trial.ang.sampleRate,2^8,(2^8)*.875,[],[],[],[1,30]);

%figure,subplot(211),bar(linspace(1,4,64),histc(max(log10(ya(Trial.stc{'g'},fa>5&fa<14,2,2)),[],2),linspace(1,4,64)),'histc')
%subplot(212),bar(linspace(1,4,64),histc(max(log10(ya(Trial.stc{'l'},fa>5&fa<14,2,2)),[],2),linspace(1,4,64)),'histc')

%figure,subplot(211),bar(linspace(1,4,128),histc(log10(Trial.xyz(Trial.stc{'g'},7,3)),linspace(1,4,128)),'histc')
%       subplot(212),bar(linspace(1,4,128),histc(log10(Trial.xyz(Trial.stc{'l'},7,3)),linspace(1,4,128)),'histc')

states = 'wgl';
nsts = numel(states);
nchan = numel(chans);
 figure,
for s = 1:nsts
x = [lbang(Trial.stc{states(s)},:)];
%[Co,f] = Comodugram(x,2^11,Trial.lfp.sampleRate,[1,120],2^10);
[Co,f] = Comodugram(x,2^9,Trial.ang.sampleRate,[1,30],2^8);

for i =1:nchan
subplot2(nchan,nsts,i,s),
imagesc(f,f,Co(:,:,i,nchan+1)'),axis xy,
if i==1,title(Trial.stc{states(s)}.label),end
if i==1&s==1,ylabel([ 'Channel: ',num2str(chans(i))]),end
caxis([-.5,.5])
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height

end

end




%% Figure 10 - Behavior Label Confusion Matrix



%xlim([1.731139e+05,2.254813e+05])
%xlim([2.862903e+05,3.185484e+05])
%xlim([6.137097e+05,6.330645e+05])
%disp(sprintf('xlim([%s,%s])',xlim));











%% Figure 11 - RHM distributions



MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%MTAConfiguration('/data/data/gravio','absolute');

%sname = 'jg05-20120315';
%sname = 'jg05-20120309';
%sname = 'jg04-20120129';
%sname = 'jg04-20120130';
%sname = 'co01-20140222';
%sname = 'jg05-20120317';

sname = 'jg05-20120310';
Trial = MTATrial(sname,'all');
%Trial.ang.load(Trial);
Trial.xyz.load(Trial);
Trial.xyz.filter(gausswin(9)./sum(gausswin(9)));

rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = Trial.com(rb);
Trial.addMarker(Trial.xyz,'hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
Trial.addMarker(Trial.xyz,'fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(gausswin(61)./sum(gausswin(61)),hcom),[1,3,2]));

ang = Trial.ang.copy;
ang.create(Trial);
ang.data = [ang(:,4,5,3),ang(:,7,11,3)];
%ang.resample(Trial.lfp);
ang = ang.data;
ang(isnan(ang(:,1)),:)=repmat(nanmean(ang),[sum(isnan(ang(:,1))),1]);

wang = WhitenSignal(ang,[],1);

[ya,fa,ta,pha,fsta] = mtchglong(wang,2^9,Trial.ang.sampleRate,2^8, ...
                                2^8*.875,[],[],[],[1,30]);



chans = [71,81,86,96];
Trial.lfp.load(Trial,chans);
wlfp = WhitenSignal(Trial.lfp.data,[],1);
yl=[];
for i = 1:Trial.lfp.size(2),
    [yl(:,:,i),fl,tl] = mtchglong(wlfp(:,i),2^12,Trial.lfp.sampleRate,2^11,2^11*0.875,[],[],[],[1,40]);
end



figure,
sp = [];
for i = 1:numel(chans)
sp(i) = subplot(6,1,i);
imagesc(tl+2^11/Trial.lfp.sampleRate,fl,log10(yl(:,:,i)'));axis xy,caxis([1,3])
end
sp(6) = subplot(6,1,5);
imagesc(ta+2^8/Trial.ang.sampleRate,fa,log10(ya(:,:,1,1)'));axis xy,caxis([-5,-2])
sp(7) = subplot(6,1,6);
imagesc(ta+2^8/Trial.ang.sampleRate,fa,log10(ya(:,:,2,2)'));axis xy,caxis([-5,-2])
linkaxes(sp,'xy');



figure,imagesc(tl,fl,log10(yl(:,:,1)'));axis xy





%% Figure - 12 Time Shifted Comodugram rhm and lfp


MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%MTAConfiguration('/data/data/gravio','absolute');

%sname = 'jg05-20120315';
sname = 'jg05-20120310';
%sname = 'jg05-20120309';
%sname = 'jg04-20120129';
%sname = 'jg04-20120130';
%sname = 'co01-20140222';

%sname = 'jg05-20120317';
chans = 68:4:95;

Trial = MTATrial(sname,'all');
%Trial.ang.load(Trial);
Trial.xyz.load(Trial);
Trial.lfp.load(Trial,chans);
Trial.lfp.resample(Trial.xyz);

rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = Trial.com(rb);
Trial.addMarker(Trial.xyz,'hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
Trial.addMarker(Trial.xyz,'fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(gausswin(61)./sum(gausswin(61)),hcom),[1,3,2]));

ang = Trial.ang.copy;
ang.create(Trial);
ang.data = ang(:,7,11,3);
ang = ang.data;
ang(isnan(ang))=nanmean(ang);
bang = ang;
lbang = WhitenSignal([Trial.lfp.data,bang],[],1);
%lbang = MTADlfp([],[],lbang,Trial.ang.sampleRate);


states = [-400,-250,-175,-70,-20,0,20,70,175,250,400];
nsts = numel(states);
nchan = numel(chans);
 figure,
for s = 1:nsts
lb = MTADlfp([],[],lbang(:,1:end-1),Trial.ang.sampleRate);
ab = MTADlfp([],[],circshift(lbang(:,end),states(s)),Trial.ang.sampleRate);
[Co,f] = Comodugram([lb(Trial.stc{'l'},:),ab(Trial.stc{'l'},:)],2^9,Trial.ang.sampleRate,[1,30],2^8);

for i =1:nchan
subplot2(nchan,nsts,i,s),
imagesc(f,f,Co(:,:,i,nchan+1)'),axis xy,
if i==1,title(Trial.stc{'l'}.label),end
if i==1&s==1,ylabel([ 'Channel: ',num2str(chans(i))]),end
caxis([-.65,.65])
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height

end

end


%% Figure 13 - Marker to marker distances
MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%MTAConfiguration('/data/data/gravio','absolute');
sname = 'jg05-20120310';
Trial = MTATrial(sname,'all');
Trial.load('ang');
Trial.load('xyz');

figure
for m = 1:8,
for o = 1:8,
subplot2(8,8,m,o)
hist(Trial.ang(Trial.ang(:,m,o,3)~=0,m,o,3),500)
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
end
end


%% Figure 14 - Spectral power for each behavior
%sname = 'jg05-20120317';
sname = 'jg05-20120310';
%sname = 'jg05-20120309';
disp = 0;
%chans = [1:2:8];
chans = [65:1:96];
marker = 'spine_lower'

Trial = MTATrial(sname,'all');
Trial.ang.load(Trial);
Trial.xyz.load(Trial);
Trial.lfp.load(Trial,chans);
Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;

wlfp = WhitenSignal(Trial.lfp.data,[],1);

matlabpool open 12

tl=[];
fl=[];
yl=[];
spectral.nfft = 2^11;
spectral.window = 2^10;
parfor i = 1:Trial.lfp.size(2),
    [yl(:,:,i),fl(:,i),tl(:,i)] = mtchglong(wlfp(:,i),spectral.nfft,Trial.lfp.sampleRate,spectral.window,spectral.window*0.875,[],[],[],[1,40]);
end
fl = fl(:,1);
tl = tl(:,1);
yld = MTADlfp('data',yl,'sampleRate',1/diff(tl(1:2,1)));
yld.data(yld.data==0)=nan;
yld.data = log10(yld.data);
yld.data = (yld.data-repmat(nanmedian(yld.data),[yld.size(1),1,1]))./repmat(nanstd(yld.data),[yld.size(1),1,1]);
tshift = round(spectral.window/2/Trial.lfp.sampleRate*yld.sampleRate);




sts='r'
evt=1;
figure,cnt=1;
for i=linspace(-2,2,40).*yld.sampleRate
imagesc(fl,1:13,sq(nanmean(yld(round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)+i-tshift),:,:),1))'),
caxis([-1,1])
text( 35,13,num2str((i)./yld.sampleRate),'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7])
pause(.2)
frms(cnt) = getframe;
cnt=cnt+1;
end

prm.fps = 5;
prm.loop = inf;
makeGif(frms,'/gpfs01/sirota/bach/homes/gravio/figures/spect_rear_on_theta.gif',prm);


th=[];
fh=[];
yh=[];
spectral.nfft = 2^9;
spectral.window = 2^7;
spectral.freq = [40,120];
parfor i = 1:Trial.lfp.size(2),
    [yh(:,:,i),fh(:,i),th(:,i)] = mtchglong(wlfp(:,i),spectral.nfft,Trial.lfp.sampleRate,spectral.window,spectral.window*0.875,[],[],[],spectral.freq);
end
fh = fh(:,1);
th = th(:,1);
yhd = MTADlfp('data',yh,'sampleRate',1/diff(th(1:2,1)));
yhd.data(yhd.data==0)=nan;
yhd.data = log10(yhd.data);
yhd.data = (yhd.data-repmat(nanmedian(yhd.data),[yhd.size(1),1,1]))./repmat(nanstd(yhd.data),[yhd.size(1),1,1]);
tshift = round(spectral.window/2/Trial.lfp.sampleRate*yhd.sampleRate);

sts='r'
evt=2;
figure,cnt=1;
for i=linspace(-2,2,200).*yhd.sampleRate
imagesc(fh,1:yhd.size(3),sq(nanmean(yhd(round(Trial.stc{sts,yhd.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)+i-tshift),:,:),1))'),
caxis([-.51,.51])
text( 80,30,num2str(i./yhd.sampleRate),'HorizontalAlignment','center','BackgroundColor',[.7 .9 .7])
pause(.2)
frms(cnt) = getframe;
cnt=cnt+1;
end

prm.fps = 10;
prm.loop = inf;
makeGif(frms,'/gpfs01/sirota/bach/homes/gravio/figures/spect_rear_off_gamma.gif',prm);


figure,
s = 'r'
sp1 = subplot(121);hold on
sp2 = subplot(122);hold on
colors= jet(numel(chans))';

for c = colors
plot(sp1,fl,-.5*find(ismember(colors',c','rows'))+nanmean(log10(yld(Trial.stc{s,yld.sampleRate}.data(2:end-1,:),:,ismember(colors',c','rows')))),'color',c)
plot(sp2,fh,-.21*find(ismember(colors',c','rows'))+nanmean(log10(yhd(Trial.stc{s,yhd.sampleRate},:,ismember(colors',c','rows')))),'color',c)
end


figure,
%subplot(121),
chan = 15;
st1 = 'l'; st2 = 'g'; 
st1 = 'r'; st2 = 'w';
boundedline(fl,mean(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),2*std(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),'-r',fl,mean(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),2*std(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,chan)),'-b','alpha')

figure,
%subplot(121),
st1 = 'l'; st2 = 'g';
boundedline(fl,mean(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),2*std(yld(Trial.stc{st1,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),'-r',fl,mean(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),2*std(yld(Trial.stc{st2,yld.sampleRate}.data(2:end-1,:)-tshift,:,3)),'-b','alpha')


subplot(122),

boundedline(fh,mean(log10(yhd(Trial.stc{'r',yhd.sampleRate}.data(2:end-1,:),:,3))),std(log10(yhd(Trial.stc{'r',yhd.sampleRate}.data(2:end-1,:),:,3))),'-r',fh,mean(log10(yhd(Trial.stc{'w',yhd.sampleRate}.data(2:end-1,:),:,3))),std(log10(yhd(Trial.stc{'w',yhd.sampleRate}.data(2:end-1,:),:,3))),'-b','alpha')


evt = 1;
sts = 'r';
figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yld.data(:,fl<12&fl>6,:),2)),round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate),round(4*yld.sampleRate),0),2))')

evt = 2;
sts = 'r';
figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yld.data(:,fl<25&fl>18,:),2)),round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate),round(4*yld.sampleRate),0),2))')


evt = 1;
sts = 'r';
figure,imagesc(sq(nanmean(GetSegs(sq(nanmean(yhd.data(:,fh<120&fh>80,:),2)),round(Trial.stc{sts,yhd.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yhd.sampleRate),round(4*yhd.sampleRate),0),2))')

tdr = sq(mean(yl(:,fl>6&fl<12,:),2)./(mean(yl(:,fl<5,:),2)+mean(yl(:,fl<18&fl>13,:),2))/2);
tdr = MTADlfp('data',tdr,'sampleRate',1/diff(tl(1:2,1)));


tdr.data = (tdr.data-repmat(nanmedian(tdr.data),[tdr.size(1),1,1]))./repmat(nanstd(tdr.data),[tdr.size(1),1,1]);


evt = 1;
sts = 'r';
figure,imagesc(sq(nanmean(GetSegs(tdr.data,round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-2*yld.sampleRate-tshift/yld.sampleRate),round(4*yld.sampleRate),0),2))')

sts = 'w';
[U,V,D] = svd(cov(yld(round(Trial.stc{sts,yld.sampleRate}.data(diff(Trial.stc{sts}.data,1,2)>200,evt)-tshift/yld.sampleRate),:,21)));
figure,imagesc(fl,fl,U),axis xy


%% Figure # Unit state place field formation

sname = 'jg05-20120317';
Trial = MTATrial(sname,'all');
Trial.stc.updateMode('auto_wbhr');
Trial.stc.load;
states = 'rwgl';

for i = 1:numel(states),
    Trial.spk.create(Trial,Trial.xyz.sampleRate,states(i));
    spk{i} = Trial.spk.copy;
    sts{i} = Trial.stc{states(i)}.copy;
    sts{i}.cast('TimeSeries');
endn

clrs = 'rbcm';

unit = 71;
figure
hold on
for i = 1:numel(states),
Lines(spk{i}(unit),[],clrs(i));
end
for i = 1:numel(states),
plot(sts{i}.data*2+i/4,clrs(i))
end
ylim([-3,4])

plot([sqrt(sum(sq([Trial.xyz(:,7,1)-75,Trial.xyz(:,7,2)+50]).^2,2))<100]-2)

%plot((Trial.xyz(:,7,3)-Trial.xyz(:,1,3))./max(Trial.xyz(:,7,3)-Trial.xyz(:,1,3)),'r')
