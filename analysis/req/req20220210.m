;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.


configure_default_args();
MjgER2016_load_data();

sampleRate = 250;
rbins = 6;
tbins = 16;
phzEdgs = linspace(0,2*pi,tbins+1);
phzCntr = mean([phzEdgs(1:end-1);phzEdgs(2:end)]);
acsdl = repmat({zeros([0,26])},[rbins,tbins]);

for tind = [17,18,19,21,22,23];
Trial = Trials{tind};
meta = sessionList(tind);
stc = copy(Trial.stc);

xyz = preproc_xyz(Trial,'trb',sampleRate);
vxy = vel(filter(copy(xyz),'ButFilter',4,2.4,'low'),{'spine_lower','nose','hcom'},[1,2]);
rhm = fet_rhm(Trial,sampleRate);
rhmPhz = phase(rhm,[5,12]);

lfp = Trial.load('lfp',[65:96]);
phz = load_theta_phase(Trial,                          ...
                       rhm,                            ...
                       meta.subject.channelGroup.theta,...
                       meta.subject.correction.thetaPhase);

% $$$ pch = fet_HB_pitchB(Trial,sampleRate);
% $$$ pch.data = pch.data(:,1);
% $$$ pchn = copy(pch);
% $$$ pchn.filter('ButFilter',4,[3,14],'bandpass');
% $$$ pchn.data = circshift(pchn(:,1),-1)-circshift(pchn(:,1),1);
% $$$ pchn.filter('ButFilter',4,[3,14],'bandpass');
% $$$ pchn.data = circshift(pchn(:,1),-1)-circshift(pchn(:,1),1);
% $$$ pchnPhz = phase(pchn,[5,12]);
% $$$ rhmPhz = phase(pchn,[5,12]);

shift = 0;
shift = align_depth_probe_csd(Trial,meta,lfp);

channelInds = [2:27]+shift;

flfp = filter(copy(lfp),'ButFilter',4,[0.8,35],'bandpass');



csdl = copy(lfp);
csdl.data = csdl.data(:,1:29);
for c = 1:29;
    %for c = 1:32;
    csdl.data(:,c) = (flfp(:,c)-flfp(:,c+3))./100;
    %csdl.data(:,c) = circshift(flfp(:,c),-1)-circshift(flfp(:,c),1);
end


rcsdl = resample(copy(csdl),rhm);

% $$$ tlfp = copy(lfp);
% $$$ tlfp.data = tlfp.data(:,meta.subject.channelGroup.theta-64);
% $$$ tlfp.resample(rhm);

phzInds = discretize(phz.data,phzEdgs);
rhmPhzInds = discretize(rhmPhz.data,linspace(-pi,pi,rbins+1));
velInds = discretize(vxy(:,2),logspace(0,2,rbins+1));
sind = logical(get(resample(cast([stc{states{5}}],'TimeSeries'),rhm),'data')) & vxy(:,2)<12 & vxy(:,2)>2;
% $$$ 
% $$$ mThetaLfpProf = [];
% $$$ for pt = 1:tbins,
% $$$     ind = phzInds==pt & sind;
% $$$     mThetaLfpProf(pt) = mean(tlfp(ind));
% $$$ end
% $$$ figure,plot([mThetaLfpProf,mThetaLfpProf])


for pr = 1:rbins,
    for pt = 1:tbins,
        ind = phzInds==pt & rhmPhzInds==pr & sind;
        ind(ind) = 1==(rem(1:sum(ind),2));
        %ind(ind) = randn([sum(ind),1])>0;
        acsdl{pr,pt} = cat(1,acsdl{pr,pt},rcsdl(ind,channelInds));
    end
end

end

phzTstat = [];
phzDf = [];
R = 3;
L = 6;
for pt = 1:tbins    
    nR = size(acsdl{R,pt},1);
    nL = size(acsdl{L,pt},1);
phzTstat(pt,:) = (mean(acsdl{R,pt}-mean(acsdl{L,pt})))...
                     ./(sqrt(((nR-1).*std(acsdl{R,pt}).^2+(nL-1).*std(acsdl{L,pt}).^2)./(nL+nR-2))...
                     .*sqrt(1/nL+1/nR));
        phzDf(pt) = nL+nR-2;
end

phzPval = [];
for c = 1:size(phzTstat,2)
    for pt = 1:tbins
        phzPval(pt,c) = 1-tcdf(abs(phzTstat(pt,c)),phzDf(pt));
    end
end
figure,imagesc(phzPval'); caxis([0,0.001]);
colormap('jet');


macsdl = repmat({[]},[tbins,1]);
for pt = 1:tbins
    for pr = 1:rbins
    macsdl{pt} = cat(1,macsdl{pt},acsdl{pr,pt});
end
end

mRandAcsdl = [];
for pt = 1:tbins
    for pr = 1:rbins
    mRandAcsdl(:,pt,pr) = median(macsdl{pt}(randsample(size(macsdl{pt},1),2000),:));
    end
end

dacsdl = [];

mOriAcsdl = [];
for pt = 1:tbins
    for pr = 1:rbins
        mOriAcsdl(:,pt,pr) = median(acsdl{pr,pt});
    end
end


sOriAcsdl = [];
for pt = 1:tbins
    for pr = 1:rbins
        sOriAcsdl(:,pt,pr) = std(acsdl{pr,pt});
    end
end

for iter = 1:100
    mRandAcsdl = [];
    for pt = 1:tbins
        for pr = 1:rbins
            mRandAcsdl(:,pt,pr) = median(macsdl{pt}(randsample(size(macsdl{pt},1),2000),:));
        end
    end
    dacsdl(:,:,:,iter) = mOriAcsdl-mRandAcsdl;
end


[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','portrait',[],1.5,1.5,0.05,0.05);
titleText = 'Theta and RHM';
vindLabels = {'All Speeds','Speeds <10cm/s','Speeds >10cm/s'}
rhmPhzLabels = {'[-\pi , 2-\pi/3]','[-2\pi/3 , -\pi/3]','[-\pi/3 , 0]',...
                '[0 , \pi/3]','[\pi/3 , 2*\pi/3]','[2*\pi/3 , \pi]'};


for pr = 1:rbins
    sax = subplot2(5,rbins,1,pr);
        imagesc(mOriAcsdl(:,:,pr));
        colormap(sax,'jet');
        caxis([-1000,500])
        %caxis([-26,26])
        if pr == rbins
            cax = colorbar(sax);
            cax.Position(1) = sum(sax.Position([1,3]))+0.01;
        end
        title({'RHM phase ',rhmPhzLabels{pr}})
    sax = subplot2(5,rbins,2,pr);
        imagesc(log10(sOriAcsdl(:,:,pr)));
        colormap(sax,'jet');
        caxis([1,3])
        %caxis([0,1.5]);
        if pr == rbins
            cax = colorbar(sax);
            cax.Position(1) = sum(sax.Position([1,3]))+0.01;
        end
    
    sax = subplot2(5,6,3,pr);
        imagesc(mean(dacsdl(:,:,pr,:),4));
        colormap(sax,'jet');
        %caxis([-5,5])
        caxis([-150,150])
        if pr == rbins
            cax = colorbar(sax);
            cax.Position(1) = sum(sax.Position([1,3]))+0.01;
        end
        
    sax = subplot2(5,6,4,pr);
        imagesc(std(dacsdl(:,:,pr,:),[],4));
        colormap(sax,'jet');    
        %caxis([0,0.75])
        caxis([0,10])
        if pr == rbins
            cax = colorbar(sax);
            cax.Position(1) = sum(sax.Position([1,3]))+0.01;
        end
        
    sax = subplot2(5,6,5,pr);    
        phzTstat = [];
        phzDf = [];
        R = 6;
        L = pr;
        for pt = 1:tbins    
            nR = size(acsdl{R,pt},1);
            nL = size(acsdl{L,pt},1);
        phzTstat(pt,:) = (mean(acsdl{R,pt}-mean(acsdl{L,pt})))...
                             ./(sqrt(((nR-1).*std(acsdl{R,pt}).^2+(nL-1).*std(acsdl{L,pt}).^2)./(nL+nR-2))...
                             .*sqrt(1/nL+1/nR));
                phzDf(pt) = nL+nR-2;
        end
        phzPval = [];
        for c = 1:size(phzTstat,2)
            for pt = 1:tbins
                phzPval(pt,c) = 1-tcdf(abs(phzTstat(pt,c)),phzDf(pt));
            end
        end
        imagesc(phzPval'); 
        colormap(sax,'bone');
        caxis([0,0.0001]);
        if pr == rbins
            cax = colorbar(sax);
            cax.Position(1) = sum(sax.Position([1,3]))+0.01;
            sax.YColor = [1,0,0];
            sax.XColor = [1,0,0];
        end
    
end


