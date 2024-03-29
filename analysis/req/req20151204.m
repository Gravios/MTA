%tsne mapping of spike on to 2D plane
cd /storage/gravio/data/project/vr_exp/Ed10-20140820/

% Load spike ClusterIds
[aClu,nclu] = LoadClu('Ed10-20140820.clu.3');
cind = ismember(aClu,1:40);
aClu = aClu(ismember(aClu,1:40));

% Load spike Waveforms and then exclude those not in the target clusters
aSpk = permute(LoadSpk('Ed10-20140820.spk.3',8,52),[3,1,2]);
aSpk = aSpk(cind,:,:);

% Load spike Features and then exclude those not in the target clusters
[aFet,nfet] = LoadFet('Ed10-20140820.fet.3');
aFet = aFet(cind,:);

% Select subset of spike for t-SNE
skp = 60;
Fet = aFet(1:skp:end,mod(1:25,3)~=0&1:25~=25); 
Clu = aClu(1:skp:end);
Spk = aSpk(1:skp:end,:,:);

% Create Color map for spike groups 
c = jet(40);
msc = c(Clu,:);

% Map Feature space to 2d space using t-SNE
figure(1203923)
mx = tsne(Fet,msc,2,8,100);

%figure(1203923)
figure(1)
xlm = xlim(gca);
ylm = ylim(gca);

mclu = 12;
hfig = figure(2); clf; 
% Plot Clusters (red=noise,blue=clusters,cyan=selected cluster)
subplot(121);hold on;
plot(mx(Clu~=1,1),mx(Clu~=1,2),'.b');
plot(mx(Clu==1,1),mx(Clu==1,2),'.r');
plot(mx(Clu==mclu,1),mx(Clu==mclu,2),'.c');
xlim(xlm);
ylim(ylm);
% Plot average spike Waveform of the seleceted cluster
subplot(122)
plot(bsxfun(@plus,sq(mean(Spk(Clu==mclu,:,:)))',2000:2000:16000))


pos = map_feature_to_tsne_space(Fet,tsneFet,tsneMap)
comptSneMap = zeros([numel(aClu),2]);
tic
for s = 1:numel(aClu),
    [~,mind] = sort(sqrt(sum(bsxfun(@minus,Fet,aFet(s,mod(1:25,3)&1:25~=25)).^2,2)));
    comptSneMap(s,:) = mean(mx(mind(1:4),:));
    if ~mod(s,10000),toc,disp([num2str(s) ' of ' num2str(numel(aClu))]),tic,end
end
toc

plot(bsxfun(@plus,sq(Spk(22,:,:))',2000:2000:16000))
plot(bsxfun(@plus,sq(Spk(mind(3),:,:))',2000:2000:16000))




for mclu = 2:36,
    hfig = figure(2); clf; 
    % Plot Clusters (red=noise,blue=clusters,cyan=selected cluster)
    subplot(121);hold on;
    plot(comptSneMap(aClu==1,1),comptSneMap(aClu==1,2),'.r');
    plot(comptSneMap(aClu~=1,1),comptSneMap(aClu~=1,2),'.b');
    plot(comptSneMap(aClu==mclu,1),comptSneMap(aClu==mclu,2),'.c');
    xlim(xlm);
    ylim(ylm);
    % Plot average spike Waveform of the seleceted cluster
    subplot(122)
    plot(bsxfun(@plus,sq(mean(aSpk(aClu==mclu,:,:)))',2000:2000:16000))
    drawnow;
    saveas(hfig,'/storage/gravio/figures/req20151204','png')
    pause(.1);
end



















% $$$ 
% $$$ function cluster_viewer(Clu,Spk,Fet)
% $$$ 
% $$$     
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ function update_cluster_display(hfig,clu)
% $$$ mclu = 9;
% $$$ hfig = figure(2); clf; 
% $$$ % Plot Clusters (red=noise,blue=clusters,cyan=selected cluster)
% $$$ subplot(121);hold on;
% $$$ plot(mx(Clu~=1,1),mx(Clu~=1,2),'.b');
% $$$ plot(mx(Clu==1,1),mx(Clu==1,2),'.r');
% $$$ plot(mx(Clu==mclu,1),mx(Clu==mclu,2),'.c');
% $$$ xlim(xlm);
% $$$ ylim(ylm);
% $$$ % Plot average spike Waveform of the seleceted cluster
% $$$ subplot(122)
% $$$ plot(bsxfun(@plus,sq(mean(Spk(Clu==mclu,:,:)))',2000:2000:16000))
% $$$     
% $$$ end

