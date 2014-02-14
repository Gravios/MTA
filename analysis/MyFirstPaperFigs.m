MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
 
Trial = MTATrial('jg05-20120317','all');
Trial.xyz.load(Trial);
Trial.filter('xyz',gausswin(31)./sum(gausswin(31)));


%% CRAP Figure 1 - Comparison between head and body motion

dmax = 100;dmin =0.1;
dxys = Trial.vel(1:Trial.xyz.size(2),[1,2]);
for i = 1:Trial.xyz.size(2),
    dxysc = clip(dxys(:,i)*10,dmin,dmax);
    dxysc(isnan(dxysc))=dmin;
    [dxyscu(:,i) dxyscux(:,i)] = MakeUniformDistr(dxysc);
end


dxyscu = MTADxyz([],[],dxyscu,Trial.xyz.sampleRate,[],[],Trial.xyz.model);
dxyscux = MTADxyz([],[],dxyscux,Trial.xyz.sampleRate,[],[],Trial.xyz.model);    

tdxyscu = dxyscu.copy;
tdxyscux = dxyscux.copy;
tdxyscu.data = dxyscu(Trial.stc{'t'},:);
tdxyscux.data = dxyscux(Trial.stc{'t'},:);

marpairs = {{{'spine_lower'},{'head_front'}},{{'pelvis_root'},{'head_back'}},{{'spine_upper'},{'head_back'}},{{'head_front'},{'head_back'}}};
nbins = 100;
nticks = 10;
for i = 1:size(marpairs,2),
    figure
    hist2([tdxyscu(:,marpairs{i}{1}) tdxyscu(:,marpairs{i}{2})],100,100)
    caxis([0 100])
    mdxyscux1 = reshape(tdxyscux(1+mod(size(tdxyscux,1),100):end,Trial.xyz.model.gmi(marpairs{i}{1})),[],nbins);
    mdxyscux2 = reshape(tdxyscux(1+mod(size(tdxyscux,1),100):end,Trial.xyz.model.gmi(marpairs{i}{2})),[],nbins);
    for j = 1:nticks
        xticks{j} = sprintf('%2.2f',mdxyscux1(round(size(mdxyscux1,1)/2),j*nticks));
        yticks{j} = sprintf('%2.2f',mdxyscux2(round(size(mdxyscux2,1)/2),j*nticks));
    end
    xl = xlim;
    yl = ylim;
    set(gca,'XTick',linspace(xl(1),xl(2),10),'XTickLabel',xticks)
    set(gca,'YTick',linspace(yl(1),yl(2),10),'YTickLabel',yticks)
end

set(gca,'XTickLabel',sprintf('%2.2f\n', ...
                             mdxyscux1(round(size(mdxyscux1,1)/2),nticks:nticks:nbins)))
set(gca,'YTickLabel',sprintf('%2.2f\n', ...
                             mdxyscux2(round(size(mdxyscux2,1)/2),nticks:nticks:nbins)))
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

MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%MTAConfiguration('/data/data/gravio','absolute');
sname = 'jg05-20120317';
sname = 'jg05-20120310';
sname = 'jg05-20120309';
Trial = MTATrial(sname,'all');
Trial.xyz.load(Trial);

Trial.load('nq');
Trial.lfp.load(Trial,[71:3:84]);



tbp = ButFilter(Trial.lfp.data(:,:),3,[6,12]./(Trial.lfp.sampleRate/2),'bandpass');
tbp_hilbert = Shilbert(tbp);
tbp_phase = phase(tbp_hilbert);
tbp_phase = MTADlfp([],[],tbp_phase,Trial.lfp.sampleRate,Trial.lfp.sync.copy,Trial.lfp.origin);


state = 'hwalk';

tbp_phase.resample(Trial.xyz);
Trial.spk.create(Trial,Trial.xyz.sampleRate,state);

%pfs = MTAAknnpfs(Trial,[],state,0,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20]);

%pfs = MTAAknnpfs(Trial,[],state,0,'numIter',1,'ufrShufBlockSize',0,'binDims',[20,20],'distThresh',100,'nNearestNeighbors',120);
pfs = MTAApfs(Trial,[],state,1,[],[20,20],[1.3,1.3]);

units = pfs.data.clu(:)';

% Get the expected ufr for each xy 
wpmr = zeros(Trial.xyz.size(1),numel(units));
wpmr = ones(Trial.xyz.size(1),numel(units));
[~,indx] = min(abs(repmat(pfs.adata.bins{1}',Trial.xyz.size(1),1)-repmat(Trial.xyz(:,7,1),1,numel(pfs.adata.bins{1}))),[],2);
[~,indy] = min(abs(repmat(pfs.adata.bins{2}',Trial.xyz.size(1),1)-repmat(Trial.xyz(:,7,2),1,numel(pfs.adata.bins{2}))),[],2);
indrm = sub2ind(pfs.adata.binSizes',indx,indy);


for unit = units,
    rateMap = reshape(pfs.data.rateMap(:,unit),pfs.adata.binSizes');
    %rateMap = pfs.plot(unit);
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

figure,
for unit = units
clf
res = Trial.spk(unit);
rind = SplitIntoBursts(res,3);
res = res(unique(rind));

if isempty(res),continue,end
if numel(res)>20000||numel(res)<200,continue,end

drzspk = DRZ(res,unit);
phzspk = tbp_phase(res,1);
gind = ~isnan(drzspk)&~isnan(phzspk);


subplot2(2,2,[1,2],1),pfs.plot(unit);
hold on,plot(pmp(unit,1),pmp(unit==units,2),'w*')

subplot2(2,2,1,2),plot(drzspk(gind),circ_rad2ang(phzspk(gind)),'.')
hold on,     plot(drzspk(gind),circ_rad2ang(phzspk(gind))+360,'.')

subplot2(2,2,2,2)
hist2([[drzspk(gind);drzspk(gind)],[circ_rad2ang(phzspk(gind));circ_rad2ang(phzspk(gind))+360]],20,15)



pause(.1)
waitforbuttonpress
end
 
 



%% End Figure 4
