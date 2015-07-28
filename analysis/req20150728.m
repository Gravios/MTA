
Trial = MTATrial('jg05-20120317')

xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);


% FXYZ Filtered Marker Positions {Low Pass 2.5 Hz}
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);



%%% Video start

rind = Trial.stc{'r'}(11,:);
ind = rind + [-60,60];
ind = ind(1):ind(2);

xlm = [min(min(xyz(ind,{'spine_lower','head_top'},1))),...
       max(max(xyz(ind,{'spine_lower','head_top'},1)))]+[-20,20];
ylm = [min(min(xyz(ind,{'spine_lower','head_top'},2))),...
       max(max(xyz(ind,{'spine_lower','head_top'},2)))]+[-20,20];
zlm = [min(min(xyz(ind,{'spine_lower','head_top'},3))),...
       max(max(xyz(ind,{'spine_lower','head_top'},3)))]+[-20,20];


hax = [];
hfig = figure(384883);clf
hfig.Position = [212 142 1372 743];

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


rstart = line(fet(rind(1),1),fet(rind(1),2));
rstart.Marker = '^';
rstart.MarkerFaceColor = 'g';
rstart.MarkerEdgeColor = 'g';
rstart.MarkerSize = 12;
rstart.Visible = 'off';

rtraj = animatedline;
rtraj.Marker = '.';
rtraj.MarkerEdgeColor = 'm';
rtraj.MarkerFaceColor = 'm';
rtraj.MarkerSize = 12;
rtraj.addpoints(fet(ind(1),1),fet(ind(1),2));

rstop = line(fet(rind(2),1),fet(rind(2),2));
rstop.Marker = 'v';
rstop.MarkerFaceColor = 'r';
rstop.MarkerEdgeColor = 'r';
rstop.MarkerSize = 12;
rstop.Visible = 'off';


legend({'rearing on set','rearing trajectory','rearing off set'})


Rex = VideoWriter(fullfile(Trial.spath,[Trial.filebase,'-rearing_ex.avi']),'Uncompressed AVI');
%Rex.Quality = 100;
%Rex.Colormap = colormap;
Rex.open;

for i = ind,
    cla(hax(1))
    rtraj.addpoints(fet(i,1),fet(i,2));
    plotSkeleton(Trial,xyz,i,'surface',ang,'hax',hax(1));
    frm = getframe(hfig);
    Rex.writeVideo(frm.cdata);
    %pause(.1)
    if i == rind(1), rstart.Visible = 'on'; end
    if i == rind(2), rstop.Visible = 'on'; end
end


Rex.close;


%%% Video End


ind = Trial.stc{'a-r'};
hist2(fet(ind,:),0:.02:(pi/2),-.04:.001:.04);
caxis([0,200]);
