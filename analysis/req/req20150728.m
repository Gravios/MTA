
set(0,'defaultAxesFontSize',14,...
      'defaultTextFontSize',14)
%encoding = 'Motion JPEG AVI';
encoding = 'Uncompressed AVI';
tag = '_raw';

% req20150728
% create video of rear 

Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);


% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

Trial.load('stc','hand_labeled_rev3_jg');

%%% Video start

% START First example Rear




hax = [];
hfig = figure(384883);clf
hfig.Position = [212 142 1372 743];
rstart = {};rtraj = {};rstop = {};



rind = Trial.stc{'r'}(11,:);
ind = rind + [-60,60];
ind = ind(1):ind(2);

xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-20,20];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-20,20];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];



hax(1) = subplot2(8,8,[1:8],[1:4]);
hold on,
plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
xlim(hax(1),xlm)
ylim(hax(1),ylm)
zlim(hax(1),zlm)
%view(hax(1),-34,10)
view(hax(1),-150,12)


hax(2) = subplot2(8,8,[2:7],[6:8]);cla
fet = Trial.xyz.copy;
fet.data = ([fang(:,1,4,2),[diff(fang(:,3,4,2));0]*fang.sampleRate]);
aind = Trial.stc{'a'};
hist2(fet(aind,:),0:.02:(pi/2),[-.04:.001:.04].*fang.sampleRate);
caxis([0,200]);
xlabel('BLBU_p_i_t_c_h (rad)');
ylabel('d(BLBU_p_i_t_c_h)/dt (rad/s)');





rtraj{1} = animatedline;
rtraj{1}.Parent = hax(2);
rtraj{1}.Marker = '.';
rtraj{1}.MarkerEdgeColor = 'm';
rtraj{1}.MarkerFaceColor = 'm';
rtraj{1}.MarkerSize = 12;
rtraj{1}.addpoints(fet(ind(1),1),fet(ind(1),2));

rstart{1} = line(fet(rind(1),1),fet(rind(1),2));
rstart{1}.Parent = hax(2);
rstart{1}.Marker = '^';
rstart{1}.MarkerFaceColor = 'g';
rstart{1}.MarkerEdgeColor = 'g';
rstart{1}.MarkerSize = 12;
rstart{1}.Visible = 'off';

rstop{1} = line(fet(rind(2),1),fet(rind(2),2));
rstop{1}.Parent = hax(2);
rstop{1}.Marker = 'v';
rstop{1}.MarkerFaceColor = 'r';
rstop{1}.MarkerEdgeColor = 'r';
rstop{1}.MarkerSize = 12;
rstop{1}.Visible = 'off';


legend({'rearing trajectory','rearing onset','rearing offset'})


% OPEN Video File
Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1',tag,'.avi']),encoding);
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1.mj2']),'Archival');
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1',tag,'.avi']),'Uncompressed AVI');
Rex.open;

% DRAW Trajectory through phase space
for i = ind,
    cla(hax(1))
    rtraj{1}.addpoints(fet(i,1),fet(i,2));
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
    if i == rind(1), rstart{1}.Visible = 'on'; end
    if i == rind(2), rstop{1}.Visible = 'on'; end
    drawnow
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end


% PAUSE for a fraction of a second 
for i = 1:10,
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end

rstart{1}.MarkerSize = 6;

rtraj{1}.MarkerSize = 6;
rtraj{1}.MarkerEdgeColor = [.6,.6,.6];
rtraj{1}.MarkerFaceColor = [.6,.6,.6];

rstop{1}.MarkerSize = 6;


% END First example


% START Second example Rear
rind = Trial.stc{'r'}(53,:);
ind = rind + [-60,60];
ind = ind(1):ind(2);

% SET new boundaries for the skeleton axes 
xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-40,40];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-40,40];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];

cla(hax(1))
plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
xlim(hax(1),xlm)
ylim(hax(1),ylm)
zlim(hax(1),zlm)
view(hax(1),-86,18)

% SET JPDF trajectories
rtraj{2} = animatedline;
rtraj{2}.Parent = hax(2);
rtraj{2}.Marker = '.';
rtraj{2}.MarkerEdgeColor = 'm';
rtraj{2}.MarkerFaceColor = 'm';
rtraj{2}.MarkerSize = 12;
rtraj{2}.addpoints(fet(ind(1),1),fet(ind(1),2));

rstart{2} = line(fet(rind(1),1),fet(rind(1),2));
rstart{2}.Parent = hax(2);
rstart{2}.Marker = '^';
rstart{2}.MarkerFaceColor = 'g';
rstart{2}.MarkerEdgeColor = 'g';
rstart{2}.MarkerSize = 12;
rstart{2}.Visible = 'off';

rstop{2} = line(fet(rind(2),1),fet(rind(2),2));
rstop{2}.Parent = hax(2);
rstop{2}.Marker = 'v';
rstop{2}.MarkerFaceColor = 'r';
rstop{2}.MarkerEdgeColor = 'r';
rstop{2}.MarkerSize = 12;
rstop{2}.Visible = 'off';


for i = ind,
    cla(hax(1))
    rtraj{2}.addpoints(fet(i,1),fet(i,2));
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
    if i == rind(1), rstart{2}.Visible = 'on'; end
    if i == rind(2), rstop{2}.Visible = 'on'; end
    drawnow
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end


% PAUSE for a fraction of a second 
for i = 1:10,
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end

rstart{2}.MarkerSize = 6;

rtraj{2}.MarkerSize = 6;
rtraj{2}.MarkerEdgeColor = [.6,.6,.6];
rtraj{2}.MarkerFaceColor = [.6,.6,.6];

rstop{2}.MarkerSize = 6;


% FILL trajectories of other rears


for s = numel(rstart)+1:10,
    rind = Trial.stc{'r'}(s,:);
    ind = rind + [-60,60];
    ind = ind(1):10:ind(2);


    rtraj{s} = animatedline;
    rtraj{s}.Parent = hax(2);
    rtraj{s}.Marker = '.';
    rtraj{s}.MarkerEdgeColor = 'm';
    rtraj{s}.MarkerFaceColor = 'm';
    rtraj{s}.MarkerSize = 12;
    rtraj{s}.addpoints(fet(ind(1),1),fet(ind(1),2));

    rstart{s} = line(fet(rind(1),1),fet(rind(1),2));
    rstart{s}.Parent = hax(2);
    rstart{s}.Marker = '^';
    rstart{s}.MarkerFaceColor = 'g';
    rstart{s}.MarkerEdgeColor = 'g';
    rstart{s}.MarkerSize = 12;
    rstart{s}.Visible = 'off';

    rstop{s} = line(fet(rind(2),1),fet(rind(2),2));
    rstop{s}.Parent = hax(2);
    rstop{s}.Marker = 'v';
    rstop{s}.MarkerFaceColor = 'r';
    rstop{s}.MarkerEdgeColor = 'r';
    rstop{s}.MarkerSize = 12;
    rstop{s}.Visible = 'off';


    for i = ind,
        rtraj{s}.addpoints(fet([i,i+10],1),fet([i,i+10],2));
        if i >= rind(1), rstart{s}.Visible = 'on'; end
        if i >= rind(2), rstop{s}.Visible = 'on'; end
        drawnow
        frm = getframe(hfig);
        Rex.writeVideo(frm.cdata);
    end

    rstart{s}.MarkerSize = 6;
    
    rtraj{s}.MarkerSize = 6;
    rtraj{s}.MarkerEdgeColor = [.6,.6,.6];
    rtraj{s}.MarkerFaceColor = [.6,.6,.6];

    rstop{s}.MarkerSize = 6;

end

% PAUSE for a fraction of a second 
for i = 1:10,
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end



Rex.close;


%%% Video End



% DECOMPOSE JPDF

mfig = figure(943882);
subplot2(2,6,1,[3,4))
aind = Trial.stc{'a'};
hist2(fet(aind,:),0:.02:(pi/2),[-.04:.001:.04].*fang.sampleRate);
caxis([0,200]);
xlabel('BLBU_{pitch} (rad)');
ylabel('d(BLBU_{pitch})/dt (rad/s)');
colorbar
for c = fliplr(20:5:200),caxis([0,c]),pause(.1),end
    
figure,
aind = Trial.stc{'a-r'};
hist2(fet(aind,:),0:.02:(pi/2),[-.04:.001:.04].*fang.sampleRate);
caxis([0,20]);
xlabel('BLBU_{pitch} (rad)');
ylabel('d(BLBU_{pitch})/dt (rad/s)');

figure,
aind = Trial.stc{'r'};
hist2(fet(aind,:),0:.02:(pi/2),[-.04:.001:.04].*fang.sampleRate);
caxis([0,20]);
xlabel('BLBU_{pitch} (rad)');
ylabel('d(BLBU_{pitch})/dt (rad/s)');




%% VIDEO Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REARING  pBMBU vs d(pBMBU)/dt phase space
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');
fang = create(MTADang,Trial,fxyz);
Trial.load('stc','hand_labeled_rev3_jg');


OwnDir = '/storage/gravio/ownCloud/SimonsFoundationTalk/';
FigDir = 'rearing';
mkdir(fullfile(OwnDir,FigDir));
FigName = [Trial.filebase,'-rearing_ex1_BMBU_singleExample'];

% START First example Rear

hax = [];
hfig = figure(384883);clf
hfig.Position = [212 142 1372 743];
rstart = {};rtraj = {};rstop = {};


rind = Trial.stc{'r',fxyz.sampleRate}(11,:);
ind = rind + [-40,40];
ind = ind(1):2:ind(2);

xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-20,20];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-20,20];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];



hax(1) = subplot2(8,8,[1:8],[1:4]);
hold on,
plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
xlim(hax(1),xlm)
ylim(hax(1),ylm)
zlim(hax(1),zlm)
%view(hax(1),-34,10)
view(hax(1),-150,12)


edx = linspace(-pi/4,pi/2,75);
edy = linspace([[-.04,.04].*fang.sampleRate,75]);
hax(2) = subplot2(8,8,[2:7],[6:8]);cla
fet = Trial.xyz.copy;
fet.data = ([fang(:,3,4,2),[diff(fang(:,3,4,2));0]*fang.sampleRate]);
aind = Trial.stc{'a'};
hist2(fet(aind,:),edx,edy);
caxis([0,200]);
xlabel('Upper Spine Pitch (rad)');
ylabel('d(Upper Spine Pitch)/dt (rad/s)');





rtraj{1} = animatedline;
rtraj{1}.Parent = hax(2);
rtraj{1}.Marker = '.';
rtraj{1}.MarkerEdgeColor = 'm';
rtraj{1}.MarkerFaceColor = 'm';
rtraj{1}.MarkerSize = 12;
rtraj{1}.addpoints(fet(ind(1),1),fet(ind(1),2));

rstart{1} = line(fet(rind(1),1),fet(rind(1),2));
rstart{1}.Parent = hax(2);
rstart{1}.Marker = '^';
rstart{1}.MarkerFaceColor = 'g';
rstart{1}.MarkerEdgeColor = 'g';
rstart{1}.MarkerSize = 12;
rstart{1}.Visible = 'off';

rstop{1} = line(fet(rind(2),1),fet(rind(2),2));
rstop{1}.Parent = hax(2);
rstop{1}.Marker = 'v';
rstop{1}.MarkerFaceColor = 'r';
rstop{1}.MarkerEdgeColor = 'r';
rstop{1}.MarkerSize = 12;
rstop{1}.Visible = 'off';


legend({'rearing trajectory','rearing onset','rearing offset'})


% OPEN Video File
Rex = VideoWriter(fullfile(OwnDir,FigDir,[FigName,'',tag,'.avi']),encoding);
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1.mj2']),'Archival');
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1',tag,'.avi']),'Uncompressed AVI');
Rex.open;

% DRAW Trajectory through phase space
for i = ind,
    cla(hax(1))
    rtraj{1}.addpoints(fet(i,1),fet(i,2));
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
    if i >= rind(1), rstart{1}.Visible = 'on'; end
    if i >= rind(2), rstop{1}.Visible = 'on'; end
    drawnow
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end


% PAUSE for a fraction of a second 
for i = 1:10,
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end

rstart{1}.MarkerSize = 6;

rtraj{1}.MarkerSize = 6;
rtraj{1}.MarkerEdgeColor = [.6,.6,.6];
rtraj{1}.MarkerFaceColor = [.6,.6,.6];

rstop{1}.MarkerSize = 6;


% END First example


% START Second example Rear
% $$$ rind = Trial.stc{'r'}(53,:);
% $$$ ind = rind + [-40,40];
% $$$ ind = ind(1):ind(2);
% $$$ 
% $$$ % SET new boundaries for the skeleton axes 
% $$$ xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
% $$$        max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-40,40];
% $$$ ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
% $$$        max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-40,40];
% $$$ zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
% $$$        max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];
% $$$ 
% $$$ cla(hax(1))
% $$$ plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
% $$$ xlim(hax(1),xlm)
% $$$ ylim(hax(1),ylm)
% $$$ zlim(hax(1),zlm)
% $$$ view(hax(1),-86,18)
% $$$ 
% $$$ % SET JPDF trajectories
% $$$ rtraj{2} = animatedline;
% $$$ rtraj{2}.Parent = hax(2);
% $$$ rtraj{2}.Marker = '.';
% $$$ rtraj{2}.MarkerEdgeColor = 'm';
% $$$ rtraj{2}.MarkerFaceColor = 'm';
% $$$ rtraj{2}.MarkerSize = 12;
% $$$ rtraj{2}.addpoints(fet(ind(1),1),fet(ind(1),2));
% $$$ 
% $$$ rstart{2} = line(fet(rind(1),1),fet(rind(1),2));
% $$$ rstart{2}.Parent = hax(2);
% $$$ rstart{2}.Marker = '^';
% $$$ rstart{2}.MarkerFaceColor = 'g';
% $$$ rstart{2}.MarkerEdgeColor = 'g';
% $$$ rstart{2}.MarkerSize = 12;
% $$$ rstart{2}.Visible = 'off';
% $$$ 
% $$$ rstop{2} = line(fet(rind(2),1),fet(rind(2),2));
% $$$ rstop{2}.Parent = hax(2);
% $$$ rstop{2}.Marker = 'v';
% $$$ rstop{2}.MarkerFaceColor = 'r';
% $$$ rstop{2}.MarkerEdgeColor = 'r';
% $$$ rstop{2}.MarkerSize = 12;
% $$$ rstop{2}.Visible = 'off';
% $$$ 
% $$$ 
% $$$ for i = ind,
% $$$     cla(hax(1))
% $$$     rtraj{2}.addpoints(fet(i,1),fet(i,2));
% $$$     plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
% $$$     if i >= rind(1), rstart{2}.Visible = 'on'; end
% $$$     if i >= rind(2), rstop{2}.Visible = 'on'; end
% $$$     drawnow
% $$$     frm = getframe(hfig);
% $$$     Rex.writeVideo(frm.cdata);
% $$$ end
% $$$ 
% $$$ 
% $$$ % PAUSE for a fraction of a second 
% $$$ for i = 1:10,
% $$$     frm = getframe(hfig);
% $$$     Rex.writeVideo(frm.cdata);
% $$$ end
% $$$ 
% $$$ rstart{2}.MarkerSize = 6;
% $$$ 
% $$$ rtraj{2}.MarkerSize = 6;
% $$$ rtraj{2}.MarkerEdgeColor = [.6,.6,.6];
% $$$ rtraj{2}.MarkerFaceColor = [.6,.6,.6];
% $$$ 
% $$$ rstop{2}.MarkerSize = 6;
% $$$ % END First example


% FILL trajectories of other rears


for s = numel(rstart)+1:20,
    rind = Trial.stc{'r'}(s,:);
    ind = rind + [-40,40];
    ind = ind(1):10:ind(2);


    rtraj{s} = animatedline;
    rtraj{s}.Parent = hax(2);
    rtraj{s}.Marker = '.';
    rtraj{s}.MarkerEdgeColor = 'm';
    rtraj{s}.MarkerFaceColor = 'm';
    rtraj{s}.MarkerSize = 12;
    rtraj{s}.addpoints(fet(ind(1),1),fet(ind(1),2));

    rstart{s} = line(fet(rind(1),1),fet(rind(1),2));
    rstart{s}.Parent = hax(2);
    rstart{s}.Marker = '^';
    rstart{s}.MarkerFaceColor = 'g';
    rstart{s}.MarkerEdgeColor = 'g';
    rstart{s}.MarkerSize = 12;
    rstart{s}.Visible = 'off';

    rstop{s} = line(fet(rind(2),1),fet(rind(2),2));
    rstop{s}.Parent = hax(2);
    rstop{s}.Marker = 'v';
    rstop{s}.MarkerFaceColor = 'r';
    rstop{s}.MarkerEdgeColor = 'r';
    rstop{s}.MarkerSize = 12;
    rstop{s}.Visible = 'off';


    for i = ind,
        rtraj{s}.addpoints(fet([i,i+10],1),fet([i,i+10],2));
        if i >= rind(1), rstart{s}.Visible = 'on'; end
        if i >= rind(2), rstop{s}.Visible = 'on'; end
        drawnow
        frm = getframe(hfig);
        Rex.writeVideo(frm.cdata);
    end

    rstart{s}.MarkerSize = 6;
    
    rtraj{s}.MarkerSize = 6;
    rtraj{s}.MarkerEdgeColor = [.6,.6,.6];
    rtraj{s}.MarkerFaceColor = [.6,.6,.6];

    rstop{s}.MarkerSize = 6;

end

% PAUSE for a fraction of a second 
for i = 1:10,
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end



Rex.close;




hfig.PaperPositionMode = 'auto';
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));


%%% Video End


%% VIDEO Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TURNING  pBMBU vs d(pBMBU)/dt phase space

% LOAD Trial Data
Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1.5,'low');
fang = create(MTADang,Trial,fxyz);
Trial.load('stc','hand_labeled_rev3_jg');


OwnDir = '/storage/gravio/ownCloud/SimonsFoundationTalk/';
FigDir = 'Turning';
mkdir(fullfile(OwnDir,FigDir));
FigName = [Trial.filebase,'-Turning_ex1_BLBU_singleExample'];

% START First example Turn
hax = [];
hfig = figure(384883);clf
hfig.Position = [10 10 1372 743];
rstart = {};rtraj = {};rstop = {};


rind = Trial.stc{'n',fxyz.sampleRate}(11,:);
ind = rind + [-40,40];
ind = ind(1):2:ind(2);

xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-20,20];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-20,20];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];



hax(1) = subplot2(8,8,[1:8],[1:4]);
hold on,
plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
xlim(hax(1),xlm)
ylim(hax(1),ylm)
zlim(hax(1),zlm)
%view(hax(1),-34,10)
view(hax(1),-150,12)


edx = linspace(-pi/4,pi/2,75);
edy = linspace([[-.04,.04].*fang.sampleRate,75]);
hax(2) = subplot2(8,8,[2:7],[6:8]);cla
fet = Trial.xyz.copy;
fet.data = [circ_dist(circ_mean(atan2(circshift(fxyz(:,:,1),-20)-circshift(fxyz(:,:,1),20),...
                               circshift(fxyz(:,:,2),-20)-circshift(fxyz(:,:,2),20)),[],2),...
                      fang(:,1,4,1)),...
            circ_dist(circshift(fang(:,1,4,1),-20),circshift(fang(:,1,4,1),20))];
aind = Trial.stc{'a'};
hist2(fet(aind,:),edx,edy);
caxis([0,200]);
xlabel('Upper Spine Pitch (rad)');
ylabel('d(Upper Spine Pitch)/dt (rad/s)');

cda = circ_dist(circ_mean(atan2(circshift(fxyz(:,:,1),-20)-circshift(fxyz(:,:,1),20),circshift(fxyz(:,:,2),-20)-circshift(fxyz(:,:,2),20)),[],2),fang(:,1,4,1));



rtraj{1} = animatedline;
rtraj{1}.Parent = hax(2);
rtraj{1}.Marker = '.';
rtraj{1}.MarkerEdgeColor = 'm';
rtraj{1}.MarkerFaceColor = 'm';
rtraj{1}.MarkerSize = 12;
rtraj{1}.addpoints(fet(ind(1),1),fet(ind(1),2));

rstart{1} = line(fet(rind(1),1),fet(rind(1),2));
rstart{1}.Parent = hax(2);
rstart{1}.Marker = '^';
rstart{1}.MarkerFaceColor = 'g';
rstart{1}.MarkerEdgeColor = 'g';
rstart{1}.MarkerSize = 12;
rstart{1}.Visible = 'off';

rstop{1} = line(fet(rind(2),1),fet(rind(2),2));
rstop{1}.Parent = hax(2);
rstop{1}.Marker = 'v';
rstop{1}.MarkerFaceColor = 'r';
rstop{1}.MarkerEdgeColor = 'r';
rstop{1}.MarkerSize = 12;
rstop{1}.Visible = 'off';


legend({'rearing trajectory','rearing onset','rearing offset'})


% OPEN Video File
Rex = VideoWriter(fullfile(OwnDir,FigDir,[FigName,'',tag,'.avi']),encoding);
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1.mj2']),'Archival');
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1',tag,'.avi']),'Uncompressed AVI');
Rex.open;

% DRAW Trajectory through phase space
for i = ind,
    cla(hax(1))
    rtraj{1}.addpoints(fet(i,1),fet(i,2));
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
    if i >= rind(1), rstart{1}.Visible = 'on'; end
    if i >= rind(2), rstop{1}.Visible = 'on'; end
    drawnow
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end


% PAUSE for a fraction of a second 
for i = 1:10,
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end

rstart{1}.MarkerSize = 6;

rtraj{1}.MarkerSize = 6;
rtraj{1}.MarkerEdgeColor = [.6,.6,.6];
rtraj{1}.MarkerFaceColor = [.6,.6,.6];

rstop{1}.MarkerSize = 6;


% END First example


% START Second example Turn
% $$$ rind = Trial.stc{'r'}(53,:);
% $$$ ind = rind + [-40,40];
% $$$ ind = ind(1):ind(2);
% $$$ 
% $$$ % SET new boundaries for the skeleton axes 
% $$$ xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
% $$$        max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-40,40];
% $$$ ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
% $$$        max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-40,40];
% $$$ zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
% $$$        max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];
% $$$ 
% $$$ cla(hax(1))
% $$$ plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
% $$$ xlim(hax(1),xlm)
% $$$ ylim(hax(1),ylm)
% $$$ zlim(hax(1),zlm)
% $$$ view(hax(1),-86,18)
% $$$ 
% $$$ % SET JPDF trajectories
% $$$ rtraj{2} = animatedline;
% $$$ rtraj{2}.Parent = hax(2);
% $$$ rtraj{2}.Marker = '.';
% $$$ rtraj{2}.MarkerEdgeColor = 'm';
% $$$ rtraj{2}.MarkerFaceColor = 'm';
% $$$ rtraj{2}.MarkerSize = 12;
% $$$ rtraj{2}.addpoints(fet(ind(1),1),fet(ind(1),2));
% $$$ 
% $$$ rstart{2} = line(fet(rind(1),1),fet(rind(1),2));
% $$$ rstart{2}.Parent = hax(2);
% $$$ rstart{2}.Marker = '^';
% $$$ rstart{2}.MarkerFaceColor = 'g';
% $$$ rstart{2}.MarkerEdgeColor = 'g';
% $$$ rstart{2}.MarkerSize = 12;
% $$$ rstart{2}.Visible = 'off';
% $$$ 
% $$$ rstop{2} = line(fet(rind(2),1),fet(rind(2),2));
% $$$ rstop{2}.Parent = hax(2);
% $$$ rstop{2}.Marker = 'v';
% $$$ rstop{2}.MarkerFaceColor = 'r';
% $$$ rstop{2}.MarkerEdgeColor = 'r';
% $$$ rstop{2}.MarkerSize = 12;
% $$$ rstop{2}.Visible = 'off';
% $$$ 
% $$$ 
% $$$ for i = ind,
% $$$     cla(hax(1))
% $$$     rtraj{2}.addpoints(fet(i,1),fet(i,2));
% $$$     plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
% $$$     if i >= rind(1), rstart{2}.Visible = 'on'; end
% $$$     if i >= rind(2), rstop{2}.Visible = 'on'; end
% $$$     drawnow
% $$$     frm = getframe(hfig);
% $$$     Rex.writeVideo(frm.cdata);
% $$$ end
% $$$ 
% $$$ 
% $$$ % PAUSE for a fraction of a second 
% $$$ for i = 1:10,
% $$$     frm = getframe(hfig);
% $$$     Rex.writeVideo(frm.cdata);
% $$$ end
% $$$ 
% $$$ rstart{2}.MarkerSize = 6;
% $$$ 
% $$$ rtraj{2}.MarkerSize = 6;
% $$$ rtraj{2}.MarkerEdgeColor = [.6,.6,.6];
% $$$ rtraj{2}.MarkerFaceColor = [.6,.6,.6];
% $$$ 
% $$$ rstop{2}.MarkerSize = 6;
% $$$ % END First example


% FILL trajectories of other rears


for s = numel(rstart)+1:20,
    rind = Trial.stc{'r'}(s,:);
    ind = rind + [-40,40];
    ind = ind(1):10:ind(2);


    rtraj{s} = animatedline;
    rtraj{s}.Parent = hax(2);
    rtraj{s}.Marker = '.';
    rtraj{s}.MarkerEdgeColor = 'm';
    rtraj{s}.MarkerFaceColor = 'm';
    rtraj{s}.MarkerSize = 12;
    rtraj{s}.addpoints(fet(ind(1),1),fet(ind(1),2));

    rstart{s} = line(fet(rind(1),1),fet(rind(1),2));
    rstart{s}.Parent = hax(2);
    rstart{s}.Marker = '^';
    rstart{s}.MarkerFaceColor = 'g';
    rstart{s}.MarkerEdgeColor = 'g';
    rstart{s}.MarkerSize = 12;
    rstart{s}.Visible = 'off';

    rstop{s} = line(fet(rind(2),1),fet(rind(2),2));
    rstop{s}.Parent = hax(2);
    rstop{s}.Marker = 'v';
    rstop{s}.MarkerFaceColor = 'r';
    rstop{s}.MarkerEdgeColor = 'r';
    rstop{s}.MarkerSize = 12;
    rstop{s}.Visible = 'off';


    for i = ind,
        rtraj{s}.addpoints(fet([i,i+10],1),fet([i,i+10],2));
        if i >= rind(1), rstart{s}.Visible = 'on'; end
        if i >= rind(2), rstop{s}.Visible = 'on'; end
        drawnow
        frm = getframe(hfig);
        Rex.writeVideo(frm.cdata);
    end

    rstart{s}.MarkerSize = 6;
    
    rtraj{s}.MarkerSize = 6;
    rtraj{s}.MarkerEdgeColor = [.6,.6,.6];
    rtraj{s}.MarkerFaceColor = [.6,.6,.6];

    rstop{s}.MarkerSize = 6;

end

% PAUSE for a fraction of a second 
for i = 1:10,
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end



Rex.close;




hfig.PaperPositionMode = 'auto';
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));


%%% Video End



%% VIDEO Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLOWRESPIRATION 

% LOAD Trial and features
Trial = MTATrial.validate('Ed05-20140529.ont.all');
ncp = fet_ncp(Trial);
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

% SETUP Video file
OwnDir = '/storage/gravio/ownCloud/SimonsFoundationTalk/';
FigDir = 'SlowRespiration';
mkdir(fullfile(OwnDir,FigDir));
FigName = [Trial.filebase,'-SlowBreathingExample_bigFont'];
Rex = VideoWriter(fullfile(OwnDir,FigDir,[FigName,'',tag,'.avi']),encoding);
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1.mj2']),'Archival');
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1',tag,'.avi']),'Uncompressed AVI');
Rex.open;


% INITIAL plots
% $$$ figure,plot((1:ncp.size(1))/ncp.sampleRate,nunity(ncp.data))
% $$$ hold on, 
% $$$ plot((1:ncp.size(1))/ncp.sampleRate,nunity(ang(:,4,5,3)).*50)
% $$$ %plot((1:ncp.size(1))/ncp.sampleRate,[0;nunity(diff(hxyz(:,3,3)))].*25)
% $$$ plot((1:ncp.size(1))/ncp.sampleRate,nunity(hxyz(:,5,3)).*100)
% $$$ Lines(Trial.stc{'s',1}(:),[],'r');


% CREATE figure
hax = [];
hfig = figure(384884);clf
hfig.Position = [212 142 1372 743];

% SET Data and Time window
ind = round(1250*ang.sampleRate):round(1266*ang.sampleRate);
timePoints = (1:ncp.size(1))/ncp.sampleRate;
features = [nunity(ncp.data),xyz(:,3,3)];

% SET Visual boundaries for skeleton animation
xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-40,40];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-40,40];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];


% PLOT first frame of skeleton animation
% $$$ hax(1) = subplot2(8,8,[1:8],[1:4]);
% $$$ hold on,
% $$$ plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
% $$$ xlim(hax(1),xlm)
% $$$ ylim(hax(1),ylm)
% $$$ zlim(hax(1),zlm)
% $$$ %view(hax(1),-34,10)
% $$$ view(hax(1),-181,4)

% CREATE axes for timeseries
%hax(2) = subplot2(8,8,[3:4],[5:8]);cla
hax(2) = subplot2(8,8,[3:4],[1:8]);cla
xlim(timePoints(ind([1,end])))
ylim([-1,1])

%hax(3) = subplot2(8,8,[6:7],[5:8]);cla
hax(3) = subplot2(8,8,[6:7],[1:8]);cla
xlim(timePoints(ind([1,end])))
ylim([62.2,63.3])



% CREATE animated lines associated with each axis
spineMiddleHeight = animatedline('Parent',hax(2));
spineMiddleHeight.LineStyle = '-';
spineMiddleHeight.Color = [1,0,0];
spineMiddleHeight.LineWidth = 1.5;

nasalCavityPressure = animatedline('Parent',hax(3));
nasalCavityPressure.LineStyle = '-';
nasalCavityPressure.Color = [.5,0,.5];
nasalCavityPressure.LineWidth = 1.5;

%set(hax(1),'Position',[0.07,0.11,0.375515463917526,0.815])

% ADD labels
%zlabel(hax(1),'Height')
xlabel(hax(3),'Time(s)')
title(hax(2),'Nasal Pressure')
ylabel(hax(2),'AU')
title(hax(3),'Spine Marker Height')
ylabel(hax(3),'mm')
% LOOP throught frames
for i = ind(1:2:end),
    %cla(hax(1))
    %plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
    spineMiddleHeight.addpoints(timePoints(i),features(i,1));
    nasalCavityPressure.addpoints(timePoints(i),features(i,2));
    drawnow
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end



Rex.close;

hfig.PaperPositionMode = 'auto';
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));







%% VIDEO Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNIFFING 

% LOAD Trial and features
Trial = MTATrial.validate('Ed05-20140529.ont.all');
ncp = fet_ncp(Trial);
rhm = fet_rhm(Trial);
fet = rhm.copy;
fet.data = [fet.data,ncp.data];
[ys,fs,ts] = fet_spec(Trial,fet,'mtchglong',false);
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');

hrb = xyz.model.rb({'head_back','head_left','head_front','head_right','head_top'});
hxyz = xyz.copy;
hxyz.data = xyz(:,hrb.ml,:);
hxyz.model = hrb;
hxyz.data = bsxfun(@minus,hxyz.data,fxyz.com(hrb));
hrb.Connections(1)=[];
hang = create(MTADang,Trial,hxyz);


% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

% SETUP Video file
OwnDir = '/storage/gravio/ownCloud/SimonsFoundationTalk/';
FigDir = 'Sniffing';
mkdir(fullfile(OwnDir,FigDir));
FigName = [Trial.filebase,'-SniffingExample_raw'];
Rex = VideoWriter(fullfile(OwnDir,FigDir,[FigName,'',tag,'.avi']),encoding);
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1.mj2']),'Archival');
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1',tag,'.avi']),'Uncompressed AVI');
Rex.open;


% INITIAL plots
% $$$ figure,plot((1:ncp.size(1))/ncp.sampleRate,nunity(ncp.data))
% $$$ hold on, 
% $$$ plot((1:ncp.size(1))/ncp.sampleRate,nunity(rhm.data).*3)
%plot((1:ncp.size(1))/ncp.sampleRate,[0;nunity(diff(hxyz(:,3,3)))].*25)
% $$$ plot((1:ncp.size(1))/ncp.sampleRate,nunity(hxyz(:,5,3)).*100)
% $$$ Lines(Trial.stc{'s',1}(:),[],'r');


% CREATE figure
hax = [];
hfig = figure(384885);clf
hfig.Position = [10 10 1372 743];
hfig.PaperPositionMode = 'auto';

% SET Data and Time window
ind = round(675*ang.sampleRate):round(690*ang.sampleRate);
timePoints = (1:ncp.size(1))/ncp.sampleRate;
features = [nunity(ncp.data),nunity(rhm.data).*3];

% SET Visual boundaries for skeleton animation
xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-40,40];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-40,40];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-40,40];


% PLOT first frame of skeleton animation
hax(1) = subplot2(8,8,[1:8],[1:4]);
hold on,
plotSkeleton(Trial,hxyz,ind(1),'surface',hang,'hax',hax(1));
xlim(hax(1),[-50,50])
ylim(hax(1),[-50,50])
zlim(hax(1),[-40,80])
%view(hax(1),-34,10)
view(hax(1),-181,4)

% $$$ % Plot Rat Skeleton
% $$$ hax(2) = subplot2(8,8,[6:8],[1:4]);
% $$$ hold on,
% $$$ plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(2));
% $$$ xlim(hax(2),xlm)
% $$$ ylim(hax(2),ylm)
% $$$ zlim(hax(2),zlm)
% $$$ %view(hax(1),-34,10)
% $$$ view(hax(2),-181,4)

% CREATE axes for timeseries

% $$$ hax(3) = subplot2(8,8,[4:5],[5:8]);cla
% $$$ xlim(timePoints(ind([1,end])))
% $$$ ylim([-1,1])

hax(4) = subplot2(8,8,[2:5],[5:8]);cla
imagesc(ts,fs(1:100),ys(:,1:100,1,2)');
axis xy
caxis([0,1])
cax = colorbar('EastOutside');
cax.Position = cax.Position+[0.05,0,0,0];
set(hax(4),'XTickLabel', {});


hax(5) = subplot2(8,8,[6:7],[5:8]);cla
ylim([-5,5]);


% CREATE animated lines associated with each axis
rhythmicHeadMotion = animatedline('Parent',hax(5));
rhythmicHeadMotion.LineStyle = '-';
rhythmicHeadMotion.Color = [1,0,0];
rhythmicHeadMotion.LineWidth = 2;

nasalCavityPressure = animatedline('Parent',hax(5));
nasalCavityPressure.LineStyle = '-';
nasalCavityPressure.Color = [.5,0,.5];
nasalCavityPressure.LineWidth = 2;

set(hax(1),'Position',get(hax(1),'Position')+[-0.07, 0, 0, 0])


% ADD labels
% $$$ zlabel(hax(2),'Height (mm)','FontSize',14)
title(hax(4),['Coherence between Nasal Cavity Pressure VS Rhythmic Head Motion'],...
      'FontSize',14)
ylabel(hax(4),'Frequency (Hz)','FontSize',14)
xlabel(hax(5),'Time(s)','FontSize',14)
ylabel(hax(5),'AU','FontSize',14)
lax = legend('Nasal Cavity Pressure','Rhythmic Head Motion','location','SouthEastOutside');
lax.FontSize = 14;
lax.Position = lax.Position + [0+.02,-0.15,0,0];

rhythmicHeadMotion.addpoints(timePoints(ind(1)-340:ind(1)),features(ind(1)-340:ind(1),1));
nasalCavityPressure.addpoints(timePoints(ind(1)-340:ind(1)),features(ind(1)-340:ind(1),2));    


% LOOP throught frames
for i = ind(1:2:end),
    cla(hax(1))
    %cla(hax(2))
    plotSkeleton(Trial,hxyz,i,'surface',hang,'hax',hax(1));
    %plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(2));
    rhythmicHeadMotion.addpoints(timePoints(i),features(i,1));
    nasalCavityPressure.addpoints(timePoints(i),features(i,2));    
    xlim(hax(4),timePoints([i-340,i+120]))    
    xlim(hax(5),timePoints([i-340,i+120]))
    drawnow
    if ismember(find(i==ind),[11,337,571]),
        print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'_frame',num2str(i),'.png']));
    end
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end

Rex.close;

print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'_frame',num2str(i),'.png']));




%% VIDEO Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WALKING

% LOAD Trial data
Trial = MTATrial.validate('jg05-20120317.cof.all');
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,4,'low');
fang = create(MTADang,Trial,fxyz);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
fvxy = fxyz.vel([1],[1,2]);

% SETUP Video file
OwnDir = '/storage/gravio/ownCloud/SimonsFoundationTalk/';
FigDir = 'Walking';
mkdir(fullfile(OwnDir,FigDir));
FigName = [Trial.filebase,'-WalkingExample' tag];
Rex = VideoWriter(fullfile(OwnDir,FigDir,[FigName,'',tag,'.avi']),encoding);
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1.mj2']),'Archival');
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1',tag,'.avi']),'Uncompressed AVI');
Rex.open;


% INITIAL plots
% $$$ figure,plot((1:ncp.size(1))/ncp.sampleRate,nunity(ncp.data))
% $$$ hold on, 
% $$$ plot((1:ncp.size(1))/ncp.sampleRate,nunity(ang(:,4,5,3)).*50)
% $$$ %plot((1:ncp.size(1))/ncp.sampleRate,[0;nunity(diff(hxyz(:,3,3)))].*25)
% $$$ plot((1:ncp.size(1))/ncp.sampleRate,nunity(hxyz(:,5,3)).*100)
% $$$ Lines(Trial.stc{'s',1}(:),[],'r');


% CREATE figure
hax = [];
hfig = figure(384884);clf
hfig.Position = [10 10 1372 743];

% SET Data and Time window
ind = 155000:158500;
ind = 1500:3100;
ind = 73000:73500;
ind = 113250:115000;
timePoints = (1:xyz.size(1))/xyz.sampleRate;
features = [[0;diff(circ_dist(fang(:,1,2,1),fang(:,3,4,1)))].*fang.sampleRate,...
            xyz(:,1,3),...
            fvxy.data];

% SET Visual boundaries for skeleton animation
xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-40,40];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-40,40];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];


% PLOT first frame of skeleton animation
hax(1) = subplot2(9,8,[6:9],[1:4]);
hold on,
plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
xlim(hax(1),xlm)
ylim(hax(1),ylm)
zlim(hax(1),zlm)
%view(hax(1),-34,10)
view(hax(1),120,4)

hax(5) = subplot2(9,8,[1:5],[1:4]);
hold on,
plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(5));
xlim(hax(5),xlm)
ylim(hax(5),ylm)
zlim(hax(5),zlm)
view(hax(5),-180,90)

% CREATE axes for timeseries
hax(2) = subplot2(9,8,[1:2],[5:8]);cla
xlim(timePoints(ind([1,end])))
ylim([-6,6])
title('d(Spine Curvature)/dt')
ylabel('rad/s')

hax(3) = subplot2(9,8,4:5,5:8);cla
xlim(timePoints(ind([1,end])))
ylim([10,70])
ylabel('Height (mm)')
title('Lower Spine Height')


hax(4) = subplot2(9,8,7:8,5:8);cla
xlim(timePoints(ind([1,end])))
title('Lower Spine speed')
ylim([0,60])
ylabel('cm/s')

% CREATE animated lines associated with each axis

spineDiffCurve = animatedline('Parent',hax(2));
spineDiffCurve.LineStyle = '-';
spineDiffCurve.Color = [1,0,0];
spineDiffCurve.LineWidth = 2;

spineLowerHeight = animatedline('Parent',hax(3));
spineLowerHeight.LineStyle = '-';
spineLowerHeight.Color = [0.5,0,0.5];
spineLowerHeight.LineWidth = 2;

spineLowerSpeed = animatedline('Parent',hax(4));
spineLowerSpeed.LineStyle = '-';
spineLowerSpeed.Color = [0,0,0];
spineLowerSpeed.LineWidth = 2;

set(hax(1),'Position',[0.07,0.11,0.375515463917526,0.347683486238532])
set(hax(5),'Position',[0.07,0.48,0.375515463917526,0.441183486238532])

% ADD labels
zlabel(hax(1),'Height')
xlabel(hax(4),'Time(s)')

% LOOP throught frames
for i = ind(1:2:end),
    cla(hax(1))
    cla(hax(5))    
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(5));
    spineDiffCurve.addpoints(timePoints(i),features(i,1));
    spineLowerHeight.addpoints(timePoints(i),features(i,2));
    spineLowerSpeed.addpoints(timePoints(i),features(i,3));
    drawnow
    if ismember(find(i==ind),[11,337,571,971,1301]),
        print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'_frame',num2str(i),'.png']));
    end    
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end



Rex.close;

hfig.PaperPositionMode = 'auto';
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));





%% VIDEO Start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TURNING

% LOAD Trial data
Trial = MTATrial.validate('jg05-20120317.cof.all');
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,4,'low');
fang = create(MTADang,Trial,fxyz);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
fvxy = fxyz.vel([1,7],[1,2]);

% SETUP Video file
OwnDir = '/storage/gravio/ownCloud/SimonsFoundationTalk/';
FigDir = 'Turning';
mkdir(fullfile(OwnDir,FigDir));
FigName = [Trial.filebase,'-TurningExample' tag];
Rex = VideoWriter(fullfile(OwnDir,FigDir,[FigName,'',tag,'.avi']),encoding);
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1.mj2']),'Archival');
%Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1',tag,'.avi']),'Uncompressed AVI');
Rex.open;


% INITIAL plots
%figure,plot(circ_dist(fang(:,1,4,1),fang(:,5,7,1)))
%figure,plot((1:ncp.size(1))/ncp.sampleRate,nunity(ncp.data))
%hold on, 
% $$$ plot((1:ncp.size(1))/ncp.sampleRate,nunity(ang(:,4,5,3)).*50)
% $$$ %plot((1:ncp.size(1))/ncp.sampleRate,[0;nunity(diff(hxyz(:,3,3)))].*25)
% $$$ plot((1:ncp.size(1))/ncp.sampleRate,nunity(hxyz(:,5,3)).*100)
% $$$ Lines(Trial.stc{'s',1}(:),[],'r');


% CREATE figure
hax = [];
hfig = figure(384884);clf
hfig.Position = [10 10 1372 743];

% SET Data and Time window
ind = 155000:158500;
ind = 1500:3100;
ind = 73000:73500;
ind = 113250:115000;
ind = 10500:11500;
timePoints = (1:xyz.size(1))/xyz.sampleRate;
features = [circ_dist(fang(:,1,4,1),fang(:,3,4,1)),...
            [0;diff(circ_dist(circshift(fang(:,1,4,1),20),circshift(fang(:,1,4,1),-20)))].*fang.sampleRate,...
            fvxy.data];

% SET Visual boundaries for skeleton animation
xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-40,40];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-40,40];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];


% PLOT first frame of skeleton animation
hax(1) = subplot2(9,8,[6:9],[1:4]);
hold on,
plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(1));
xlim(hax(1),xlm)
ylim(hax(1),ylm)
zlim(hax(1),zlm)
%view(hax(1),-34,10)
view(hax(1),120,4)

hax(5) = subplot2(9,8,[1:5],[1:4]);
hold on,
plotSkeleton(Trial,xyz,ind(1),'surface',ang,'hax',hax(5));
xlim(hax(5),xlm)
ylim(hax(5),ylm)
zlim(hax(5),zlm)
view(hax(5),-180,90)

% CREATE axes for timeseries
hax(2) = subplot2(9,8,[1:2],[5:8]);cla
xlim(timePoints(ind([1,end])))
ylim([-pi/2,pi/2])
title('Head-Body Offset')
ylabel('rad')

hax(3) = subplot2(9,8,4:5,5:8);cla
xlim(timePoints(ind([1,end])))
ylim([-5,5])
title('Angular Body Speed')
ylabel('rad/s')


hax(4) = subplot2(9,8,7:8,5:8);cla
xlim(timePoints(ind([1,end])))
title('Head and Body Speed')
ylim([0,60])
ylabel('cm/s')

% CREATE animated lines associated with each axis

spineDiffCurve = animatedline('Parent',hax(2));
spineDiffCurve.LineStyle = '-';
spineDiffCurve.Color = [1,0,0];
spineDiffCurve.LineWidth = 2;

spineLowerHeight = animatedline('Parent',hax(3));
spineLowerHeight.LineStyle = '-';
spineLowerHeight.Color = [0.5,0,0.5];
spineLowerHeight.LineWidth = 2;

spineLowerSpeed = animatedline('Parent',hax(4));
spineLowerSpeed.LineStyle = '-';
spineLowerSpeed.Color = [0,1,0];
spineLowerSpeed.LineWidth = 2;

headLowerSpeed = animatedline('Parent',hax(4));
headLowerSpeed.LineStyle = '-';
headLowerSpeed.Color = [1,0,0];
headLowerSpeed.LineWidth = 2;

set(hax(1),'Position',[0.07,0.11,0.375515463917526,0.347683486238532])
set(hax(5),'Position',[0.07,0.48,0.375515463917526,0.441183486238532])

% ADD labels
zlabel(hax(1),'Height')
xlabel(hax(4),'Time(s)')

% LOOP throught frames
for i = ind(1:2:end),
    cla(hax(1))
    cla(hax(5))    
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(5));
    spineDiffCurve.addpoints(timePoints(i),features(i,1));
    spineLowerHeight.addpoints(timePoints(i),features(i,2));
    spineLowerSpeed.addpoints(timePoints(i),features(i,3));
    headLowerSpeed.addpoints(timePoints(i),features(i,4));    
    drawnow
    if ismember(find(i==ind),[11,337,571,971]),
        print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'_frame',num2str(i),'.png']));
    end    
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
end



Rex.close;

hfig.PaperPositionMode = 'auto';
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));



