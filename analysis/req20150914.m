%% tSNE stuff
%mkdir('/storage/gravio/figures/req/req20150914');

hostPath = '/storage/gravio/figures/req/req20150914';
msr = 15; % New sample rate


Trial = MTATrial('jg05-20120317');
Trial.stc.load(Trial,'hand_labeled_rev2');

xyz = Trial.load('xyz');


man = Trial.load('fet','lsppc');
man.filter('ButFilter',3,2,'low');
man.resample(msr);


fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');

ang = create(MTADang,Trial,xyz);
tsh = 1;
afet = Trial.xyz.copy;
bfet = Trial.xyz.copy;
afet.data = circshift(xyz(:,:,[1,2]),-tsh)-circshift(fxyz(:,:,[1,2]),tsh);
afet.data = reshape(afet.data,[],2);
aft = mat2cell(afet.data,size(afet,1),[1,1]);
[afet.data,bfet.data] = cart2pol(aft{:});
afet.data = reshape(afet.data,[],xyz.size(2));
bfet.data = reshape(bfet.data,[],xyz.size(2));
m = MTADxyz('data',circ_dist(afet(:,1),ang(:,1,4,1)),'sampleRate',Trial.xyz.sampleRate);
m.data = circ_dist(circshift(m.data,-5),circshift(m.data,5));
m.data = [diff(m.data);0];
wn = 40;
% $$$ mv = m.copy;
% $$$ mv.data = circshift(1./permute(mean(m.segs(1:m.size(1),wn,nan).^2),[2,3,4,1]),-round(wn/2));
% $$$ mv.data(isnan(mv.data))=0;
% $$$ mv.resample(msr);
bfet.resample(msr);



fvelxy = xyz.vel([],[1,2]);
fvelxy.resample(msr);
fvelxy.filter('ButFilter',3,2.4,'low');

fvelz = fxyz.vel([],[3]);
fvelz.resample(msr);
fvelz.filter('ButFilter',3,2.4,'low');


fxyz.resample(msr);
dsa = create(MTADang,Trial,fxyz);


fet = fxyz.copy;
fet.data = [fxyz(:,[1,4,7],3),...
            fvelxy(:,[1,4,7]),....
            fvelz(:,5),...
            man.data,...
            log10(bfet(:,1)+1),...
            dsa(:,'spine_middle','spine_upper',2),...
            dsa(:,'spine_upper','head_back',2),...
            dsa(:,'head_back','head_front',2),...            
            dsa(:,1,4,3).*cos(dsa(:,1,4,2)),...
            abs(circ_dist(circshift(dsa(:,3,4,2),-1),circshift(dsa(:,3,4,2),1))),...
            abs(circ_dist(circshift(dsa(:,1,4,1),-1),circshift(dsa(:,1,4,1),1))),...
            abs(circ_dist(circshift(dsa(:,3,7,1),-1),circshift(dsa(:,3,7,1),1)))];
fet.data(isinf(fet(:))) = 0;


[asmat,labels,keys] =  stc2mat(Trial.stc,fet);
asmat = MTADxyz('data',asmat,'sampleRate',fet.sampleRate);
[~,asmat.data] = max(asmat.data(:,1:12),[],2);
c = jet(numel(unique(asmat)));
csmat = asmat.copy; 
csmat.data = c(csmat.data,:);

ind = Trial.stc{'a'};

mfet = nunity(fet(ind,:));
msmat = csmat(ind,:);

start = 1;
skip = 2;
stop = 45000;
no_dims = 2;

initial_dims = 5;
perplexity = 100;
ind = start:skip:stop;
mappedX = tsne(mfet(ind,:), msmat(ind,:), no_dims, initial_dims, perplexity);

figTitle = ['tSNE-msr_' num2str(msr) '-ind_' num2str(start) '_' ...
            num2str(skip) '_' num2str(stop) '-perplexity_' ...
            num2str(perplexity) '-no_dims_' num2str(no_dims)];
hfig = figure(1);

%saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '.fig']),'fig')
saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '.eps']),'epsc')
saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '.png']),'png')


hfig = figure;
plot(mappedX(:,1),mappedX(:,2),'.');
cls = ClusterPP(hfig);

mcid = asmat(Trial.stc{'a'});
mcid = mcid(ind);

figure,bar(1:12,histc(mcid,1:12),'histc');
figure,bar(1:12,histc(mcid(cls==1),1:12),'histc');

figure
hist2(mappedX,...
      linspace(min(mappedX(:,1)),max(mappedX(:,1)),100),...
      linspace(min(mappedX(:,2)),max(mappedX(:,2)),100));






mtfet =  MTADxyz('data',mfet(ind,:),'sampleRate',fet.sampleRate);
mtpos =  MTADxyz('data',mappedX,'sampleRate',fet.sampleRate);

[RateMap,Bins,distdw,distIndw]= PlotKNNPF(Trial,mtfet,mtpos,[5,5],20,5,'xy',[],[],[-110,110;-110,110]);
rmap = reshape(RateMap,numel(Bins{1}),numel(Bins{2}),[]);


fett = {};
fetd = {};

%% Feature tags and definitions
%lower spine speed
fett(end+1) = {'Height_{BL} (z-score)'};
fetd(end+1) = {'1 Hz low pass filtered height of the lower spine maker'};

fett(end+1) = {'Height_{BU} (z-score)'};
fetd(end+1) = {'1 Hz low pass filtered height of the upper spine maker'};

fett(end+1) = {'Height_{HF} (z-score)'};            
fetd(end+1) = {'1 Hz low pass filtered height of the head front maker'};

fett(end+1) = {'XY Speed_{BL} (z-score)'};
fetd(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                'the spine lower maker']};

fett(end+1) = {'XY Speed_{BU} (z-score)'};
fetd(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                'the spine upper maker']};

fett(end+1) = {'XY Speed_{HF} (z-score)'};
fetd(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                'the head front maker']};

fett(end+1) = {'Vertical Speed(flp1Hz) of Middle Spine (z-score)'};
fetd(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                'head back marker']};

fett(end+1) = {'PPC_{traj yaw} (z-score)'};
fetd(end+1) = {['1 Hz lowpass filtered Pair-wise Phase Consisistency(PPC) of the yaw of ' ...
                'trajectories of all makers along the rostro-caudal axis']};

fett(end+1) = {'bfet (z-score)'};
fetd(end+1) = {['Magnitude of the projection of lower spine trajectory  ' ...
                'onto the vecor of lower spine to upper spine']};

fett(end+1) = {'Pitch_{BMBU} (z-score)'};
fetd(end+1) = {['Pitch of spine_middle to spine_upper relative to xy ' ...
                'plane']};

fett(end+1) = {'Pitch_{BUHB} (z-score)'};
fetd(end+1) = {['Pitch of spine_upper to head_back relative to xy ' ...
                'plane']};

fett(end+1) = {'Pitch_{HBHF} (z-score)'};
fetd(end+1) = {['Pitch of head_back to head_front relative to xy ' ...
                'plane']};

fett(end+1) = {'XY Dist_{BLBU} (z-score)'};
fetd(end+1) = {['Magnitude of the projection of the vector formed ' ...
                'by the spine_lower and spine_upper markers']};

fetd(end+1) = {'d(pitch_{BMBU})/dt (z-score'};
fetd(end+1) = {'Pitch speed of the vector from spine_middle to spine_upper'};

fetd(end+1) = {'d(yaw_{BLBU})/dt (z-score'};
fetd(end+1) = {'Pitch speed of the vector from spine_middle to spine_upper'};

fetd(end+1) = {'d(yaw_{BMHF})/dt (z-score'};
fetd(end+1) = {'Pitch speed of the vector from spine_middle to head_front'};



hfig = figure(38380);
for i = 1:fet.size(2);
subplot(132); hold on

sts = Trial.stc.list_state_attrib('label');
mc = msmat(ind,:);
for nc = 1:12,
    nind = all(bsxfun(@eq,c(nc,:),mc),2);
    h = scatter(mappedX(nind,1),mappedX(nind,2),16,mc(nind,:));
    h.MarkerFaceColor = h.CData(1,:);
end
legend(sts(1:12))
xlim([min(mappedX(:,1))-5,max(mappedX(:,1))+30]);
ylim([min(mappedX(:,2))-5,max(mappedX(:,2))+5]);
daspect([1,1,1])

subplot(133),cla
[~,hc] = imagescnan({Bins{1},Bins{2},rmap(:,:,4)},[],false,true,[0,0,0]);
axis xy
daspect([1,1,1])


subplot{133)

figTitle = ['tSNE-msr_' num2str(msr) '-ind_' num2str(start) '_' ...
            num2str(skip) '_' num2str(stop) '-perplexity_' ...
            num2str(perplexity) '-no_dims_' num2str(no_dims)];


saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '.eps']),'epsc')
saveas(hfig,fullfile(hostPath,[Trial.filebase '_' figTitle '.png']),'png')

end

