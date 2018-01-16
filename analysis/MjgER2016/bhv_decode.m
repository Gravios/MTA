addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/fwdrebayesiandecodingtools/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/tempScripts/
addpath /storage/antsiro/data/lab/homes_of_alumni/marcel/scripts/aux/

Trial = MTATrial('jg05-20120310');
Trial = MTATrial('jg05-20120311');
Trial.load('stc','msnn_ppsvd_raux');
%states = Trial.stc.list_state_attrib;
states = {'loc&theta','rear&theta','pause&theta'};

Trial.load('nq');
pft = pfs_2d_theta(Trial); 
mrt = pft.maxRate;
units = select_units(Trial,18,'pyr');
units = units(mrt(pft.data.clu(units))>0.8);

for sts = 1:numel(states),
    defargs = get_default_args('MjgER2016','MTAApfs','struct');
    defargs.units = units;
    defargs.states = states{sts};
    defargs.numIter = 1000;
    defargs.halfsample = true;
    defargs.overwrite = false;
    defargs = struct2varargin(defargs);
    pfs{sts} = MTAApfs(Trial,defargs{:});
end    

[accg,tbin] = autoccg(Trial,units,'theta');

t = 1;
mRate = max(cell2mat(cf(@(p) p.maxRate, pfs)));
%units = find(sq(mRate(1,1,:))>3);

%Circ mask
width = pfs{1}.adata.binSizes(1);
height = pfs{1}.adata.binSizes(2);
radius = round(pfs{1}.adata.binSizes(1)/2)-...
         find(pfs{1}.adata.bins{1}<-420,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
mask(mask==0)=nan;



%sunits = find(sq(mRate(:))>3&sq(mRate(:))<30);
%sunits =  sunits(~ismember(sunits,[4 16 32 63 97 98 99 100 106 107]));
% $$$ sunits = units;
% $$$ 
% $$$ autoincr = false;
% $$$ hfig = figure(849274);    
% $$$ unit = sunits(1);
% $$$ while unit~=-1,
% $$$     for s = 1:numel(states)
% $$$         subplot(3,4,s)
% $$$         hold('on')
% $$$         pf = pfs{s};
% $$$         ratemap = pf.plot(unit,'isCircular',true);
% $$$         ratemap(isnan(ratemap)) = -1;
% $$$         ratemap = ratemap.*mask;
% $$$         imagesc(pf.adata.bins{1},pf.adata.bins{2},ratemap');    
% $$$         text(pf.adata.bins{1}(1)+30,pf.adata.bins{2}(end)-50,...
% $$$              sprintf('%2.1f',max(ratemap(:))),'Color','w','FontWeight','bold','FontSize',10)
% $$$         colormap([0,0,0;parula]);
% $$$ 
% $$$ 
% $$$ % $$$             plot(peakPatchCOM(t,i,unit==units,1),...
% $$$ % $$$                  peakPatchCOM(t,i,unit==units,2),'*k');
% $$$ % $$$             xlim([-600,600]),ylim([-350,350])                    
% $$$         title([pf.session.trialName ':' pf.parameters.states,': ',num2str(unit)]);
% $$$     end   
% $$$     ForAllSubplots('colormap([0,0,0;parula])');
% $$$     ForAllSubplots(['caxis([-1,',num2str(sq(mRate(unit))),'.*1.5])']);
% $$$ 
% $$$     subplot(3,4,10)
% $$$     bar(tbin,accg(:,unit));
% $$$     xlim([min(tbin),max(tbin)]);
% $$$     title([' AutoCCG: Unit ',num2str(unit)]);
% $$$ 
% $$$     unit = figure_controls(hfig,unit,sunits,autoincr);
% $$$ end
% $$$ 


sunits = units;

ratemap = zeros([numel(pfs{1}.adata.bins{1}),...
                 numel(pfs{1}.adata.bins{2}),...
                 numel(states),...
                 numel(sunits)]);

for u = 1:numel(sunits),
    for s = 1:numel(states),
        trm = pfs{s}.plot(sunits(u),'mean');
        ratemap(:,:,s,u) = trm;
        ratemap(:,:,s,u) = ratemap(:,:,s,u).*mask;
    end
end

%ratemap(isnan(ratemap)) = 0;

xyz = Trial.load('xyz');
xyz.resample(10);
ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'theta-groom-sit',sunits,1.5,true);

E = decode_bayesian_poisson(ratemap,ufr.data');

xbins = pfs{1}.adata.bins{1};
ybins = pfs{1}.adata.bins{2};


pos = nan([size(E,4),2]);
sts = nan([size(E,4),1]);
for tind = 1:size(E,4),
    try,
        xyi = LocalMinimaN(-E(:,:,:,tind),0,100);        
        pos(tind,:) = [xbins(xyi(1)),ybins(xyi(2))];
        sts(tind,1) = xyi(3);
    end
end    



figure,
plot(sq(xyz(:,7,[1])))
hold on,
plot(pos(:,1));
plot(sts.*100);

i = 25;
sts = Trial.stc{'w',10};
ind = sts.data(i,1)-10:sts.data(i,2)+5;
figure,hold('on');
plot(sq(xyz(ind,7,[1])),sq(xyz(ind,7,[2])));
plot(sq(xyz(ind(1),7,[1])),sq(xyz(ind(1),7,[2])),'*m');
plot(pos(ind,1)+2*randn([numel(ind),1]),pos(ind,2)+2*randn([numel(ind),1]),'.');
plot(pos(ind(1),1),pos(ind(1),2),'*m');
ylim([-500,500])
xlim([-500,500])

stst = zeros([size(sts,1),7])';
inind = [1:size(sts,1)].*7-7;
stst(sts(nniz(sts))+inind(nniz(sts))') = 1;

%figure,plot(sts)
sto = stc2mat(Trial.stc,xyz,states(2:8));



shl = MTADxyz('data',double(0<sto),'sampleRate',xyz.sampleRate);
ysm = MTADxyz('data',double(stst'),'sampleRate',xyz.sampleRate); 



aper = Trial.stc{'t'};
aper.cast('TimeSeries');
aper.resample(xyz);
ind = any(shl.data,2)&any(ysm.data,2)&logical(aper.data);

tcm = confmat(shl(ind,:),ysm(ind,:)); % #DEP: netlab
confusionMatrix = round(tcm./xyz.sampleRate,2);
precision = round(diag(tcm)./sum(tcm,2),4).*100;
sensitivity = round(diag(tcm)'./sum(tcm),4).*100;
accuracy = sum(diag(tcm))/sum(tcm(:));


xp = sq(xyz(:,7,[1,2]))-pos;
[th,r] = cart2pol(xp(:,1),xp(:,2));

figure,plot(r(sum(ufr.data,2)>10))
Lines([tind],[],'k');

figure,hist(r(sum(ufr.data,2)>10),100)


figure,plot(sq(xyz(:,5,[1])))
hold on,
plot(pos(:,1).*double(sum(ufr.data,2)>100));



% Theta rate map pos decode


tratemap = zeros([numel(pfs{1}.adata.bins{1}),...
                 numel(pfs{1}.adata.bins{2}),...
                 numel(sunits)]);

for u = 1:numel(sunits),
        trm = pfs{9}.plot(sunits(u));
        tratemap(:,:,u) = trm;
        tratemap(:,:,u) = tratemap(:,:,u).*mask;
end
tE = decode_bayesian_poisson(tratemap,ufr.data');
tpos = nan([size(tE,3),2]);
for tind = 1:size(tE,3),
    try,
        xyi = LocalMinimaN(-tE(:,:,tind),0,100);        
        tpos(tind,:) = [xbins(xyi(1)),ybins(xyi(2))];
    end
end    

figure,
plot(sq(xyz(:,7,[1])))
hold on,
plot(tpos(:,1).*double(sum(ufr.data,2)>100));



txp = sq(xyz(:,7,[1,2]))-tpos;
[tht,tr] = cart2pol(txp(:,1),txp(:,2));

%figure,plot(r(sum(ufr.data,2)>10))

figure,hist(r(sum(ufr.data,2)>5),100)
figure,hist(tr(sum(ufr.data,2)>5),100)


figure,plot(sq(xyz(:,5,[1])))
hold on,
plot(tpos(:,1).*double(sum(ufr.data,2)>100));





spk = Trial.spk.copy();
spk.create(Trial,xyz.sampleRate,'theta',units,'deburst');    
[mpr,mpl] = cf(@(p) p.maxRate, pfs);

figure,
plot(mpl{1}(:,1),mpl{1}(:,2),'.');
ylim([-500,500])
xlim([-500,500])
