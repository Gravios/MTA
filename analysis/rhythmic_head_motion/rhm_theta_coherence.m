configure_default_args();
MjgER2016_load_data();

sampleRate = 250;
cshifts = [0,-6500:125:6500];
scoher = repmat({zeros(38,6,3,numel(cshifts))},[1,numel(Trials)]);


for tind = 1:numel(Trials);
    Trial = Trials{tind};
    disp(Trial.filebase);

    stc = copy(Trial.stc);
    xyz = preproc_xyz(Trial,'trb',sampleRate);
    rhm = fet_rhm(Trial,sampleRate);

    Trial.lfp.filename = [Trial.name,'.lfp'];    
    tlfp = Trial.load('lfp',sessionList(tind).subject.channelGroup.theta);
    tlfp.resample(rhm);

    nsit = ~logical(get(resample(cast([stc{'s+m'}],'TimeSeries'),rhm),'data'));

    for iter = 1:numel(cshifts)
        disp(iter)
        tic
        rhmCs = copy(rhm);
        rhmCs.data(nniz(rhmCs.data)&nsit) = circshift(rhm.data(nniz(rhmCs.data)&nsit),cshifts(iter));
        elfp = copy(tlfp);
        elfp.data = cat(2,elfp.data,rhmCs.data);
        specArgsTheta = struct('nFFT',2^9,...
                               'Fs',  elfp.sampleRate,...
                               'WinLength',2^8,...
                               'nOverlap',2^8*0.5,...
                               'NW',3,...
                               'Detrend',[],...
                               'nTapers',[],...
                               'FreqRange',[1,20]);
        [mys,mfs,mts] = fet_spec(Trial,elfp,'mtcsdglong',false,[],specArgsTheta,[],true);
        if iter==1
            mvxy = vel(resample(copy(xyz),mys),{'nose'},[1,2]);    
            vind = {};
            vind{1} = true(size(mvxy(:,1)));
            vind{2} = mvxy(:,1)<10&mvxy(:,1)>2;
            vind{3} = mvxy(:,1)>10;
        end
        for sts = 1:numel(states)
            sind = logical(get(resample(cast([stc{states{sts}}],'TimeSeries'),mys),'data'));
            for vts = 1:numel(vind)
                ind = vind{vts}&sind;
                scoher{tind}(:,sts,vts,iter) = sum(abs(mys(ind,:,1,2)))./sum(sqrt(mys(ind,:,1,1).*mys(ind,:,2,2)));
            end
        end
        toc
    end
end


figure,
subplot(311);
imagesc(mts,mfs,log10(abs(mys(:,:,1,1)))');
axis('xy');
caxis([3,5]);
colormap('jet');
subplot(312);
imagesc(mts,mfs,log10(abs(mys(:,:,2,2)))');
axis('xy');
caxis([-7,-3]);
colormap('jet');
subplot(313);
imagesc(mts,mfs,(abs(mys(:,:,1,2))./sqrt(mys(:,:,1,1).*mys(:,:,2,2)))');
axis('xy');
caxis([0,1]);
colormap('jet');
linkx();


tind = 1;
figure,
hold('on');
for sts = 1:6
plot(mfs,scoher{tind}(:,sts,1,2))
end

tind = 20;
[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','landscape',[],1.5,1.5,0.05,0.05);
titleText = 'Theta and RHM';
vindLabels = {'All Speeds','Speeds <10cm/s','Speeds >10cm/s'}
for vst = 1:numel(vind)
    for sts = 1:numel(states)
        subplot2(numel(vind),numel(states),vst,sts);
        hold('on')
        plot(mfs,mean(sq(scoher{tind}(:,sts,vts,1)),2),'b');
        plot(mfs,mean(sq(scoher{tind}(:,sts,vts,2:end)),2),'k');
        plot(mfs,mean(sq(scoher{tind}(:,sts,vts,2:end)),2)+std(sq(scoher{tind}(:,sts,vts,2:end)),[],2),'r');
        plot(mfs,mean(sq(scoher{tind}(:,sts,vts,2:end)),2)-std(sq(scoher{tind}(:,sts,vts,2:end)),[],2),'r');
        title({titleText,states{sts},['for ',vindLabels{vst}]});
        xlabel('Frequency');
        ylabel('Coherence');
        ylim([0.3,0.55]);
    end
end

tind = 20;
[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','landscape',[],1.5,1.5,0.05,0.05);
titleText = 'Theta and RHM';
vindLabels = {'All Speeds','Speeds <10cm/s','Speeds >10cm/s'}

for vst = 1:numel(vind)
    for sts = 5
        subplot2(numel(vind),numel(states),vst,sts);
        hold('on')
        plot(mfs,mean(sq(scoher{tind}(:,sts,vts,1)),2),'b');
        plot(mfs,mean(sq(scoher{tind}(:,sts,vts,2:end)),2),'k');
        plot(mfs,mean(sq(scoher{tind}(:,sts,vts,2:end)),2)+std(sq(scoher{tind}(:,sts,vts,2:end)),[],2),'r');
        plot(mfs,mean(sq(scoher{tind}(:,sts,vts,2:end)),2)-std(sq(scoher{tind}(:,sts,vts,2:end)),[],2),'r');
        title({titleText,states{sts},['for ',vindLabels{vst}]});
        xlabel('Frequency');
        ylabel('Coherence');
        ylim([0.3,0.55]);
    end
end


acoher = cat(5,scoher{:});

figure,
hold('on');
for t = 1:30
plot(mfs,acoher(:,3,1,1,t));
end


figure,
hold('on');
for t = 1:30
plot(acoher(mfs<12&mfs>5,3,1,1,t),acoher(mfs<12&mfs>5,5,1,1,t));
end
xlim([0.2,0.8]);ylim([0.2,0.8])

figure,
hold('on');
for t = 1:30
plot(acoher(mfs<12&mfs>5,3,1,100,t),acoher(mfs<12&mfs>5,5,1,100,t));
end
xlim([0.2,0.8]);ylim([0.2,0.8])


plot(mfs,mean(sq(scoher{tind}(:,sts,vts,1)),2)
