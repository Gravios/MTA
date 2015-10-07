function req20150925(Trial,varargin)
[figPath,sampleRate,sts,initial_dims,perplexity] = DefaultArgs(varargin,{'/storage/gravio/figures/',15,{'walk','rear'},5,100});

SAVEFIG = true;
osts = numel(sts)+1;

%try
%    Trial.stc.load(Trial,'auto_wbhr');
%catch
Trial = labelBhv(Trial);
%end


fet = fet_tsne(Trial,sampleRate,true);


[asmat,labels,keys] =  stc2mat(Trial.stc,fet,sts);
asmat = cat(2,asmat,~any(asmat,2));
asmat = MTADxyz('data',asmat,'sampleRate',fet.sampleRate);
[~,asmat.data] = max(asmat.data(:,[1:osts-1,end]),[],2);
c = jet(numel(unique(asmat)));
csmat = asmat.copy; 
csmat.data = c(csmat.data,:);

ind = Trial.stc{'a'};

mfet = nunity(fet(ind,:));
msmat = csmat(ind,:);

start = 1;
skip = 2;
stop = size(mfet,1);
no_dims = 2;

ind = start:skip:stop;
mappedX = tsne(mfet(ind,:), msmat(ind,:), no_dims, initial_dims, perplexity);


figTitle = ['tSNE-sampleRate_' num2str(sampleRate) '-ind_' num2str(start) '_' ...
            num2str(skip) '_' num2str(stop) '-perplexity_' ...
            num2str(perplexity) '-initial_dims_' num2str(initial_dims) ...
            '-no_dims_' num2str(no_dims)];

save(fullfile(Trial.spath,[Trial.filebase '-req20150925-' figTitle '.mat']),...
     'mappedX','mfet','msmat','sts','initial_dims','perplexity');




hfig = figure(3923923),clf,hold on
sts = Trial.stc.list_state_attrib('label');
mc = msmat(ind,:);
for nc = 1:osts,
    nind = all(bsxfun(@eq,c(nc,:),mc),2);
    h = scatter(mappedX(nind,1),mappedX(nind,2),10,mc(nind,:));
    h.MarkerFaceColor = h.CData(1,:);
end
legend([sts(1:osts-1),{'unlabeled'}])
xlim([min(mappedX(:,1))-5,max(mappedX(:,1))+30]);
ylim([min(mappedX(:,2))-5,max(mappedX(:,2))+5]);
daspect([1,1,1])


reportfig(figPath, hfig, 'tsne', 'req', false,Trial.filebase,figTitle,[],SAVEFIG,'png');


mtfet =  MTADxyz('data',mfet(ind,:),'sampleRate',fet.sampleRate);
mtpos =  MTADxyz('data',mappedX,'sampleRate',fet.sampleRate);

[RateMap,Bins,distdw,distIndw]= PlotKNNPF(Trial,mtfet,mtpos,[5,5],20,5,'xy',[],[],[-110,110;-110,110]);
rmap = reshape(RateMap,numel(Bins{1}),numel(Bins{2}),[]);




fett = {};
fetd = {};

%% Feature tags and definitions
%lower spine speed
fett(end+1) = {'Height_{BL}'};
fetd(end+1) = {'1 Hz low pass filtered height of the lower spine maker'};

fett(end+1) = {'Height_{BU}'};
fetd(end+1) = {'1 Hz low pass filtered height of the upper spine maker'};

fett(end+1) = {'Height_{HF}'};            
fetd(end+1) = {'1 Hz low pass filtered height of the head front maker'};

fett(end+1) = {'XY Speed_{BL}'};
fetd(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                'the spine lower maker']};

fett(end+1) = {'XY Speed_{BU}'};
fetd(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                'the spine upper maker']};

fett(end+1) = {'XY Speed_{HF}'};
fetd(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                'the head front maker']};

fett(end+1) = {'Vertical Speed(flp1Hz) of Middle Spine'};
fetd(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                'head back marker']};

fett(end+1) = {'PPC_{traj yaw}'};
fetd(end+1) = {['1 Hz lowpass filtered Pair-wise Phase Consisistency(PPC) of the yaw of ' ...
                'trajectories of all makers along the rostro-caudal axis']};

fett(end+1) = {'bfet'};
fetd(end+1) = {['Magnitude of the projection of lower spine trajectory  ' ...
                'onto the vecor of lower spine to upper spine']};

fett(end+1) = {'Pitch_{BMBU}'};
fetd(end+1) = {['Pitch of spine_middle to spine_upper relative to xy ' ...
                'plane']};

fett(end+1) = {'Pitch_{BUHB}'};
fetd(end+1) = {['Pitch of spine_upper to head_back relative to xy ' ...
                'plane']};

fett(end+1) = {'Pitch_{HBHF}'};
fetd(end+1) = {['Pitch of head_back to head_front relative to xy ' ...
                'plane']};

fett(end+1) = {'XY Dist_{BLBU}'};
fetd(end+1) = {['Magnitude of the projection of the vector formed ' ...
                'by the spine_lower and spine_upper markers']};

fett(end+1) = {'d(pitch_{BMBU})/dt'};
fetd(end+1) = {'Pitch speed of the vector from spine_middle to spine_upper'};

fett(end+1) = {'d(yaw_{BLBU})/dt'};
fetd(end+1) = {'Pitch speed of the vector from spine_lower to spine_upper'};

fett(end+1) = {'d(yaw_{BMHF})/dt'};
fetd(end+1) = {'Yaw speed of the vector from spine_middle to head_front'};


hfig = figure(38381);
for i = 1:size(rmap,3);
    clf
    [ha,hc] = imagescnan({Bins{1},Bins{2},rmap(:,:,i)},[],false, ...
                        true,[0,0,0]);
    ylabel(hc,'mean z-score','FontName','Courier','FontSize',12);
    axis xy
    title(fett{i},'FontName','Courier','FontSize',12)
    daspect([1,1,1])
   
    reportfig(figPath, hfig, 'tsne', 'req',false,fett{i},fetd{i},[],SAVEFIG,'png',8,8);
end
