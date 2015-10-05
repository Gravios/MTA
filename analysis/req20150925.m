% Testing various parameters for tSNE 
%

initial_dims = 5;
perplexity = 100;




%% tSNE stuff
%mkdir('/storage/gravio/figures/req/req20150925');

Trial = MTATrial('Ed05-20140529','all','ont');
Trial.stc.load(Trial,'auto_wbhr');
sts = {'walk','rear'};
osts = numel(sts)+1;

figPath = '/storage/gravio/figures/labbook/';
sampleRate = 15; % New sample rate

fet = fet_tnse(Trial,msr,true);


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


figTitle = ['tSNE-msr_' num2str(msr) '-ind_' num2str(start) '_' ...
            num2str(skip) '_' num2str(stop) '-perplexity_' ...
            num2str(perplexity) '-initial_dims_' num2str(initial_dims) ...
            '-no_dims_' num2str(no_dims)];

hax = gca;
t = all(bsxfun(@eq,hax.Children.CData,[0,1,1]),2);
hax.Children.CData(t,:) = repmat([0,1,0],[sum(t),1]);
t = all(bsxfun(@eq,hax.Children.CData,[1,1,0]),2);
hax.Children.CData(t,:) = repmat([1,0,0],[sum(t),1]);

hfig = figure(1);
hfig.PaperPosition = [0,0,6,6];
saveas(hfig,fullfile(figPath,[Trial.filebase '_' figTitle '.fig']),'fig')
saveas(hfig,fullfile(figPath,[Trial.filebase '_' figTitle '.eps']),'epsc')
saveas(hfig,fullfile(figPath,[Trial.filebase '_' figTitle '.png']),'png')





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
fetd(end+1) = {'Pitch speed of the vector from spine_middle to spine_upper'};

fett(end+1) = {'d(yaw_{BMHF})/dt'};
fetd(end+1) = {'Pitch speed of the vector from spine_middle to head_front'};



hfig = figure(38380);
hfig.Position = [23 490 1893 486];
hfig.PaperPosition = [0,0,40,10];
for i = 1:fet.size(2);
    clf
    hax = subplot(131);
    cla;
    emptyAxis(hax);
    delete(hax.Children )
    ht = text(.05,.9,...
              {['Trial:           ' Trial.filebase],...
               ['StateCollection: ' Trial.stc.mode]},...
              'Interpreter','none','FontName','Courier');
    ht.FontSize = 12;

    fdesc = [{['Feature Description:']},strchp(fetd{i},40)];    
    ht = text(.05,.75,...
              fdesc,...
              'Interpreter','none','FontName','Courier');
    ht.FontSize = 12;
    
    ht = text(.05,.65,...
              ['    States:          Tot Occ (sec): '],...
              'Interpreter','none','FontName','Courier');
    ht.FontSize = 12;

    yps = fliplr(linspace(.1,.55,osts));
    for j = 1:osts,
        sts = Trial.stc.states{j};
        ht = text(.05,yps(j),...
                  ['    ' sts.label],...
                  'Interpreter','none','FontName','Courier');
        ht.FontSize = 12;
        ht = text(.45,yps(j),...
                  [num2str(sum(diff(sts.data,1,2))./sts.sampleRate)],...
                  'Interpreter','none','FontName','Courier');
        ht.FontSize = 12;
    end

    
    subplot(132);cla
    [ha,hc] = imagescnan({Bins{1},Bins{2},rmap(:,:,i)},[],false, ...
                        true,[0,0,0]);
    ylabel(hc,'mean z-score','FontName','Courier');
    axis xy
    title(fett{i},'FontName','Courier')
    daspect([1,1,1])

   
    subplot(133),cla
   % MEAN FET MAP
    hold on
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



    saveas(hfig,fullfile(hostPath,[Trial.filebase '_featureOverLay_' num2str(i) '_' figTitle '.eps']),'epsc')
    %saveas(hfig,fullfile(hostPath,[Trial.filebase '_featureOverLay_' num2str(i) '_' figTitle '.png']),'png')

end

PlotSessionErrors