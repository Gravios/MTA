
Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);


% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);



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
Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex1.avi']),'Uncompressed AVI');
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

