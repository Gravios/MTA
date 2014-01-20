function [RateMap,Bins,MRate,SI,Spar]= PlotKNNPF(Session,ufr,pos,varargin)
[binDims,SmoothingWeights,type] = DefaultArgs(varargin,{50,[],'xy'});

Session = MTATrial('jg05-20120317');
Session.xyz.load(Session);
Session = Session.filter('xyz');


ufr = Session.ufr.create(Session,Session.xyz,[],[],0.2);


mytx = Session.xyz(Session.stc{'t'},7,[1,2]);
mytu = ufr(Session.stc{'t'},:);
spind = 1:30:size(mytx,1);
mytx = mytx(spind,:);
mytu = mytu(spind,:);


mywx = Session.xyz(Session.stc{'w'},7,[1,2]);
mywu = ufr(Session.stc{'w'},:);
spind = 1:10:size(mywx,1);
mywx = mywx(spind,:);
mywu = mywu(spind,:);

mygx = Session.xyz(Session.stc{'g'},7,[1,2]);
mygu = ufr(Session.stc{'g'},:);
spind = 1:30:size(mygx,1);
mygx = mygx(spind,:);
mygu = mygu(spind,:);

mylx = Session.xyz(Session.stc{'l'},7,[1,2]);
mylu = ufr(Session.stc{'l'},:);
spind = 1:30:size(mylx,1);
mylx = mylx(spind,:);
mylu = mylu(spind,:);


myrx = Session.xyz(Session.stc{'r'},7,[1,2]);
myru = ufr(Session.stc{'r'},:);
spind = 1:30:size(myrx,1);
myrx = myrx(spind,:);
myru = myru(spind,:);



nxbins = 50;
nybins = 50;
xbins = linspace(Session.maze.boundaries(1,1),Session.maze.boundaries(1,2),nxbins)';
ybins = linspace(Session.maze.boundaries(2,2),Session.maze.boundaries(2,1),nxbins)';
nnn = 15;
pfknnmrw = nan(nxbins,nybins,ufs.size(2));
pfknnmdw = nan(nxbins,nybins);
pfknnmrg = nan(nxbins,nybins,ufs.size(2));
pfknnmdg = nan(nxbins,nybins);
pfknnmrl = nan(nxbins,nybins,ufs.size(2));
pfknnmdl = nan(nxbins,nybins);

pfknnmrr = nan(nxbins,nybins,ufs.size(2));
pfknnmdr = nan(nxbins,nybins);


%new knnpf-start

nxbins = 50;
nybins = 50;
xbins = linspace(Session.maze.boundaries(1,1),Session.maze.boundaries(1,2),nxbins)';
ybins = linspace(Session.maze.boundaries(2,2),Session.maze.boundaries(2,1),nxbins)';
nnn = 15;
pfknnmrw = nan(nxbins,nybins,ufs.size(2));
pfknnmdw = nan(nxbins,nybins);
dthresh = 70;

xy = cell(2,1);
[xy{:}]= meshgrid(xbins,ybins);

xyi = cell(2,1)
[xyi{:}] = meshgrid(1:nxbins,1:nybins);
xyi = repmat(permute(cat(3,xyi{:}),[4,3,1,2]),size(mywx,1),1);

xy = repmat(permute(cat(3,xy{:}),[4,3,1,2]),size(mywx,1),1);
mywxt = repmat(mywx,[1,1,nxbins,nybins]);
mywut = repmat(mywu,[1,1,nxbins,nybins]);
distw = sqrt(sum((mywxt-xy).^2,2));
[distdw,distIndw ] = sort(distw,1,'ascend');  

pfknnmdw = sq(median(distdw(1:nnn,:,:,:),1));


s = size(mywut);
cm = reshape(repmat([0:s(1):prod(s)-1],s(1),1),s);     
     
smywut = mywut(repmat(distIndw,1,ufr.size(2))+cm);
pfknnmrw = permute(mean(smywut(1:70,:,:,:),1),[3,4,2,1]);
pfknnmrw = permute(mean(smywut(1:nnn,:,:,:),1),[3,4,2,1]);
%pfknnmrw = permute(mean(mywut(repmat(distIndw(1:nnn,:,:,:),1,ufr.size(2))),1),[3,4,2,1]);

pfknnmrw(repmat(pfknnmdw,[1,1,ufr.size(2)])>dthresh) = nan;
pfknnmrw(repmat(pfknnmdw,[1,1,ufr.size(2)])>40) = nan;

%new knnpf-end

for y = 1:length(ybins),
    for x = 1:length(xbins),
        distw = sqrt(sum((mywx-repmat([xbins(x),ybins(y)],size(mywx,1),1)).^2,2));
        distr = sqrt(sum((myrx-repmat([xbins(x),ybins(y)],size(myrx,1),1)).^2,2));
        distg = sqrt(sum((mygx-repmat([xbins(x),ybins(y)],size(mygx,1),1)).^2,2));
        distl = sqrt(sum((mylx-repmat([xbins(x),ybins(y)],size(mylx,1),1)).^2,2));
        [~,distIndw ] = sort(distw);        
        [~,distIndr] = sort(distr);        
        [~,distIndg ] = sort(distg);        
        [~,distIndl] = sort(distl);        
        pfknnmdw(y,x) = median(distw(distIndw(1:nnn)));
        pfknnmdr(y,x) = median(distr(distIndr(1:nnn)));
        pfknnmdg(y,x) = median(distg(distIndg(1:nnn)));
        pfknnmdl(y,x) = median(distl(distIndl(1:nnn)));
for u=1:95,
        if pfknnmdw(y,x)<70,
            %pfknnmrw(y,x) = min(mywu(distIndw(1:nnn)));
            %pfknnmrw(y,x) = max(mywu(distIndw(1:nnn)));
            pfknnmrw(y,x,u) = mean(mywu(distIndw(1:nnn),u));
        end
        if pfknnmdr(y,x)<70,
            %pfknnmrr(y,x) = min(myru(distIndr(1:nnn)));
            %pfknnmrr(y,x) = max(myru(distIndr(1:nnn)));
            pfknnmrr(y,x,u) = mean(myru(distIndr(1:nnn),u));
        end
        if pfknnmdg(y,x)<70,
            %pfknnmrw(y,x) = min(mywu(distIndw(1:nnn)));
            %pfknnmrw(y,x) = max(mywu(distIndw(1:nnn)));
            pfknnmrg(y,x,u) = mean(mygu(distIndg(1:nnn),u));
        end
        if pfknnmdl(y,x)<70,
            %pfknnmrr(y,x) = min(myru(distIndr(1:nnn)));
            %pfknnmrr(y,x) = max(myru(distIndr(1:nnn)));
            pfknnmrl(y,x,u) = mean(mylu(distIndl(1:nnn),u));
        end
end
    end
end



hfig = figure;

u = 1;
set(hfig,'Name',num2str(u));
while u~=-1,

subplot(251);
hist(mywu(:,u))
xl = xlim
xlim([-0.05,xl(2)])

subplot(252);
imagescnan({xbins,ybins,pfknnmrw(:,:,u)},[],[],0,[0,0,0]),axis xy
title(['walk u ' num2str(Session.spk.map(u,:))])
text(xbins(1)+30,ybins(1)-50,sprintf('%2.1f',max(max(pfknnmrw(:,:,u)))),'Color','w','FontWeight','bold','FontSize',10)

subplot(253);
imagescnan({xbins,ybins,pfknnmrr(:,:,u)},[],[],0,[0,0,0]),axis xy
title(['rear u ' num2str(Session.spk.map(u,:))])
text(xbins(1)+30,ybins(1)-50,sprintf('%2.1f',max(max(pfknnmrr(:,:,u)))),'Color','w','FontWeight','bold','FontSize',10)

subplot(254);
imagescnan({xbins,ybins,pfknnmrg(:,:,u)},[],[],0,[0,0,0]),axis xy
title(['hwalk u ' num2str(Session.spk.map(u,:))])
text(xbins(1)+30,ybins(1)-50,sprintf('%2.1f',max(max(pfknnmrg(:,:,u)))),'Color','w','FontWeight','bold','FontSize',10)

subplot(255);
imagescnan({xbins,ybins,pfknnmrl(:,:,u)},[],[],0,[0,0,0]),axis xy
title(['lwalk u ' num2str(Session.spk.map(u,:))])
text(xbins(1)+30,ybins(1)-50,sprintf('%2.1f',max(max(pfknnmrl(:,:,u)))),'Color','w','FontWeight','bold','FontSize',10)


subplot(256);
bar(tbin,accg(:,u)),axis tight,

subplot(257);
pfw.plot(u),

subplot(258);
pfr.plot(u),

subplot(259);
pfh.plot(u)

subplot(2,5,10);
pfl.plot(u)

%print(gcf,'-dpsc2',[Trial.filebase '.knnpf_100ms_50msol_walk_rear' num2str(u) '.eps']);
    uicontrol(hfig,'String','>','units','normalized','Position',[.97, 0, .03,  1],'Callback',@(hfig,index,indArray,callback)figure_controls_gui(hfig,u,units,'forwardButton_Callback'));
    uicontrol(hfig,'String','<','units','normalized','Position',[  0, 0, .03,  1],'Callback',@(hfig,index,indArray,callback)figure_controls_gui(hfig,u,units,'backwardButton_Callback'));
    uicontrol(hfig,'String','X','units','normalized','Position',[.45, 0, .10,.07],'Callback',@(hfig,index,indArray,callback)figure_controls_gui(hfig,u,units,'exitButton_Callback'));
    uicontrol(hfig,'String','Print','units','normalized','Position',[.55, 0, .10,.07],'Callback',@(hfig,index,indArray,callback)figure_controls_gui(hfig,u,units,'printButton_Callback'));
    while strcmp(num2str(u),get(hfig,'Name'))
        pause(0.2);
    end
    u = str2double(get(hfig,'Name'));
end
 
