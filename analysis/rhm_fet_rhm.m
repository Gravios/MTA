Trial = MTATrial('Ed05-20140528');
%Trial = MTATrial('jg05-20120317');
Trial.stc.updateMode('qda_filtf1p5');,Trial.stc.load;
xyz = Trial.xyz.copy;
xyz.load(Trial);

wang = WhitenSignal([fet_rhm(Trial),fet_ncp(Trial,'chans',[1,2])],[],1);

[ys,fs,ts] = mtcsdglong(wang,2^8,Trial.ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,30]);
szy = size(ys);
padding = zeros([round(2^6/xyz.sampleRate/diff(ts(1:2))),szy(2:end)]);
ys = MTADlfp('data',cat(1,padding,ys,padding),'sampleRate',1/diff(ts(1:2)));

% Speed of the head
xyz.filter(gtwin(1,Trial.xyz.sampleRate));
vh = MTADxyz('data',xyz.vel(1,[1,2]),'sampleRate',Trial.xyz.sampleRate);

%% Bug: zeros should remain after resample, anti-alias filter seems
%% to be affecting the edges
vh.resample(ys);



%edges_labels = mat2cell(edges,2,ones(1,size(edges,2)))';
chan_labels = {'Rhythmic Head Motion','Nasal Epithelium Signal','Nasal Cavity Pressure'};

for s = 1:numel(Trial.stc.states);

sind = Trial.stc.states{s};


vhs = log10(abs(vh(sind)));
ind = ~isinf(vhs)&~isnan(vhs);
vhs = vhs(ind);
vhlim =prctile(vhs,[5,98]);
edges = linspace(vhlim(1),vhlim(2),9);
edges = [edges(1:end-1);edges(2:end)];
edges_labels = mat2cell(10.^mean(edges),1,ones(1,size(edges,2)))';
yss = ys(sind,:,:,:);
yss = yss(ind,:,:,:);
clear ind;

vsc = [];
for i = edges,
    ind = i(1)>=vhs&vhs<i(2);% & 5>=mean(log10(yss(:,fs>13&fs>5,3,3)),2);
    for j = 1:size(yss,3),
        for k = 1:size(yss,3),
            if sum(ind)~=0,
                if k~=j,
                    vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(abs(yss(ind,:,k,j)))./mean(sqrt(yss(ind,:,k,k).*yss(ind,:,j,j)));
                else
                    vsc(find(i(1)==edges(1,:)),:,k,j) = nanmean(yss(ind,:,k,j));
                end

            else
                vsc(find(i(1)==edges(1,:)),:,:) = zeros([1,size(yss,2),size(yss,3)]);
            end

        end
    end
end


for j = 1:size(yss,3),
    for k = 1:size(yss,3),
        subplot2(size(yss,3),size(yss,3),j,k),
        if j~=k
            imagesc(1:size(edges,2),fs,vsc(:,:,j,k)'),axis xy,
            title([chan_labels{j} ' & ' chan_labels{k} ' coherence'])
        else
            imagesc(1:size(edges,2),fs,log10(vsc(:,:,j,k))'),axis xy,
            title([chan_labels{j} ' PSD Binned by Body Speed'])
        end
        set(gca,'XTickLabel',cellfun(@num2str,cellfun(@transpose,edges_labels(2:2:8),'UniformOutput',false),'UniformOutput',false)');
        xlabel('Binned Body Speed cm/s')
        ylabel('Frequency Hz')
        colorbar
    end
end

reportfig('Trial','/gpfs01/sirota/bach/homes/gravio/figures', ...
          'FileName','mean_Coherence_RHM_NCP','Comment',['State: ' ...
                    Trial.stc.states{s}.label]  )

end

yp=ys.copy;
yp.resample(Trial.xyz);
phs = Trial.lfp.phase([5,13]);
pang = LocalMinima(dbang,5);
myp = mean(log10(yp(:,fs<13&fs>5,3,3)),2);

pedges = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,;...
         0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6];

apr = [];
for i = pedges,
    ind = find(i(1)<=myp & myp<i(2) & ~isnan(myp) &~isinf(myp) & myp~=0);
    if sum(ind)~=0,
        apr(find(i(1)==pedges(1,:)),:) = hist(phs(pang(ismember(pang,ind)),2),linspace(-pi,pi,15));
        apr(find(i(1)==pedges(1,:)),:) = apr(find(i(1)==pedges(1,:)),:)./sum(apr(find(i(1)==pedges(1,:)),:));
    else
        apr(find(i(1)==pedges(1,:)),:) = zeros([1,15]);
    end
end



figure,
imagesc(linspace(-pi,pi,15),pedges(2,:),apr),axis xy
caxis([0,0.12])
title('Filtered[5-13] RHM peak NCP phase Distrb')
xlabel('Phase (radians)')
ylabel('log10 NCP 5-13hz Power')





phs = Trial.lfp.phase([5,13]);
pang = LocalMinima(dbang,5);
myp = mean(log10(yp(:,fs<13&fs>5,3,3)),2);
vhp = Trial.vel(11,[1,2]);
vhp = MTADxyz('data',[0;vhp],'sampleRate',Trial.xyz.sampleRate);
 

pedges = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5;...
         0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6];
edges = [0 1 2 4  8 12 16 20;...
         1 2 4 8 12 16 20 24];
vpr = [];
vps = [];
for i = pedges,
for j = edges,
    ind = find(i(1)<=myp & myp<i(2) & ~isnan(myp) & ~isinf(myp) & myp~=0 & j(1)>=vhp.data&vhp.data<j(2) );
    if sum(ind)~=0&&sum(ismember(pang,ind))~=0,
        vpr(find(i(1)==pedges(1,:)),find(j(1)==edges(1,:)),:) = circ_mean(phs(pang(ismember(pang,ind)),2));
        [~,vps(find(i(1)==pedges(1,:)),find(j(1)==edges(1,:)),:)] = circ_mtest(phs(pang(ismember(pang,ind)),2),vpr(find(i(1)==pedges(1,:)),find(j(1)==edges(1,:)),:));
    else
        vpr(find(i(1)==pedges(1,:)),:) = zeros([1,1,1]);
    end
end
end


figure,
imagescnan({edges(2,:),pedges(2,:),vpr},[-pi,pi],1,1,[0,0,0]),axis xy
caxis([0,0.1])
title('Filtered[5-13] RHM peak NCP phase Distrb')
xlabel('Phase (radians)')
ylabel('log10 NCP 5-13hz Power')



[ya,fa,ta,pha,fsa] = mtchglong(wang,2^8,Trial.ang.sampleRate,2^7,2^7*.875,[],[],[],[1,30]);

figure,hold on
sp = [];
for i = 1:size(ya,3),
sp(i) = subplot(3,1,i);
imagesc(ta,fa,log10(ya(:,:,i,i)'))
axis xy
end
linkaxes(sp,'xy');


[yo,fo, phi]= mtchd(wang,2^8,Trial.ang.sampleRate,2^7,2^7*.875,[],'linear',[],[1,30]);

figure,imagesc(ta,fa,log10(ya(:,:,1,3)')),axis xy





