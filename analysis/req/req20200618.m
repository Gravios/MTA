% req20200618
%    Tags: egocentric, phase precession, head-body-angle
%    Status: active
%    Type: Analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: MjgER2016: bhvfields
%    Description: egocentric coordinates relative to placefield for all placefields
%    New_methods: 


global MTA_PROJECT_PATH

MjgER2016_load_data();

states = {'theta',...
          'rear',...
          'hloc',...
          'hpause',...
          'lloc',...
          'lpause',...
          'pause',...
          'walk',...
          'turn'};

% SELECT placefields with placefield away from maze wall
units = {};
units{1}  =  [15,96,99,150,207];                           % er01-20110719 CA2/3
units{2}  =  [20,75,86,99,152,157,159];                    % er01-20110721 CA2/3

units{3}  =  [95,158,173];                                 % ER06-20130612 CA1
units{4}  =  [29,31,107,126,175,197];                      % ER06-20130613 CA1
units{5}  =  [31,35,40,48,67,81,99,119,121];               % ER06-20130614 CA1

units{6}  =  [4,7,10,18,25,35,37,38,49,60,68,104];         % Ed10-20140816 CA2/3
units{7}  =  [1,10,24,33,38,57,63,64,73,82,105,108];       % Ed10-20140817 CA2/3

units{8}  =  [];                                           % jg04-20120128 CA1
units{9}  =  [];                                           % jg04-20120129 CA1
units{10} =  [];                                           % jg04-20120130 CA1
units{11} =  [];                                           % jg04-20120131 CA1
units{12} =  [];                                           % jg04-20120201 CA1

units{13} =  [];                                           % jg04-20120210 CA3
units{14} =  [];                                           % jg04-20120211 CA3
units{15} =  [];                                           % jg04-20120212 CA3
units{16} =  [];                                           % jg04-20120213 CA3

units{17} =  [7,16,20,58,81];                              % jg05-20120309 CA1
%units{18} =  [9,11,18,33,42,49,50,54,60,74,75,78,80];      % jg05-20120310 CA1
units{18} =  [9,11,33,34,42,49,50,54,55,60,61,75,78,80,81,87];% jg05-20120310 CA1
units{19} =  [27,28,33,50,53,63,66,73,150,153];            % jg05-20120311 CA1
units{20} =  [20,21,25,35,44,61,72,79,80,103,104,110,139]; % jg05-20120312 CA1
units{21} =  [6,22,24,25,37,43,61,68,77];                  % jg05-20120315 CA1
units{22} =  [13,19,30,41,42,48,56,58,61,65];              % jg05-20120316 CA1
units{23} =  [29,35,40,50,54,63,72];                       % jg05-20120317 CA1

%units{24} =  []; er01-20110722? CA3
%units{25} =  [8,10,39,42,47,56,57,81,83,92,98,100,101,102]; % Ed10-20140815 CA2/3
%units{26} =  [4,43,47,70,82,87,140,145,155,167,175,211,230]; %ER06-20130624 CA2/3
%units{27} =  [11,17,18,24,26,28,48,49,52,53,54]; % jg05-20120323 CA2/3
%units{28} =  [8,10,19,22,29,39,40,55]; % jg05-20120324 CA2/3


tind = find(~cellfun(@isempty,units));
rot = [0,0,0,0,0,0,0,0,0,0,...
       0,0,0,0,0,0,0.17, 0.17, 0.17, 0.17,...
       0.17, 0.17, 0.17];


spkvn = load_spk_vars(sessionList,Trials,units,'ego_spk_analysis',states,rot);
% spkv is old var

% CORRECT phase
tphz = spkvn.phz;
tphz = tphz + pi*double(spkvn.map(:,1)==1|spkvn.map(:,1)==2|spkvn.map(:,1)==6|spkvn.map(:,1)==7);
tphz = tphz + pi*1/4*2*double(spkvn.map(:,1)==3|spkvn.map(:,1)==4|spkvn.map(:,1)==5);
tphz = tphz + pi*1/8*2*double(spkvn.map(:,1)==8|spkvn.map(:,1)==9|spkvn.map(:,1)==10);
tphz = tphz + pi*1/8*2*double(spkvn.map(:,1)==11|spkvn.map(:,1)==12|spkvn.map(:,1)==13|spkvn.map(:,1)==14);
tphz(tphz<0) = tphz(tphz<0) + 2*pi;
tphz(tphz>2*pi) = tphz(tphz>2*pi) - 2*pi;

% plot the 2 scatter for limited phase, with hba as color
% plot the 3 scatter for x y and phase, with hba as color

figure,
ind = ismember(spkvn.map(:,1),[3,4,5,17,18,19,20,21,22,23]) ...
      & spkvn.stc(:,1)==1 ...
      &(spkvn.stc(:,1)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9);
scatter3(spkvn.ego(ind,1),   ...
         spkvn.ego(ind,15),   ...
         tphz(ind),     ...
         10,                ...
         spkvn.hba(ind),     ...
         'filled');



chba = -(spkvn.hba+0.2*double(ismember(spkvn.map(:,1),[3,4,5])) - ...
         0.2*double(ismember(spkvn.map(:,1),[8:14])));


phzBins = 0:pi/4.5:2*pi;                  phzBinInd = discretize(tphz,phzBins);
phzBins = linspace(0,2*pi,8);             phzBinInd = discretize(tphz,phzBins);
hbaBins = [-1.2,-0.4,0.4,1.2];            hbaBinInd = discretize(chba,hbaBins);
hvlBins = [-80,-4,4,80];                  hvlBinInd = discretize(spkvn.hvl,hvlBins);
hvfBins = [-4,10,20,30,40,50,60,70,80];   hvfBinInd = discretize(spkvn.hvf,hvfBins);
bmaBins = [-pi./,-pi/2,-pi/4,0,pi]

nx = numel(phzBins)-1;
ny = numel(hvfBins)-1;
nv = numel(hvlBins)-1;


figure();
plotSkeleton(Trials{t},xyz{t},1000);
hold('on');
plot(bsxfun(@plus,xyz{t}(1000,'hcom',[1]),[0,hvec(1000,1,1)]*100),...
     bsxfun(@plus,xyz{t}(1000,'hcom',[2]),[0,hvec(1000,1,2)]*100),...
     '-');
plot(bsxfun(@plus,xyz{t}(1000,'hcom',[1]),[0,hvec(1000,2,1)]*100),...
     bsxfun(@plus,xyz{t}(1000,'hcom',[2]),[0,hvec(1000,2,2)]*100),...
     '-');
plot(0,0,'*k');

    multiprod(bsxfun(@minus,...
                         [0,0],...
                         sq(xyz{t}(1000:1010,'hcom',[1,2]))),...
                  hvec(1000:1010,:,:),2,[2,3])
             

% DST color
% place field
% theta phase, hx, hy, hba

headBodyAng = [xyz{t}(:,'spine_upper',[1,2])-xyz{t}(:,'bcom',[1,2]),...
               xyz{t}(:,'nose',[1,2])-xyz{t}(:,'hcom',[1,2])];
headBodyAng = sq(bsxfun(@rdivide,headBodyAng,sqrt(sum(headBodyAng.^2,3))));
headBodyAng = cart2pol(headBodyAng(:,:,1),headBodyAng(:,:,2));
headBodyAng = circ_dist(headBodyAng(:,2),headBodyAng(:,1));
headBodyAng = MTADfet.encapsulate(Trials{tind},-(headBodyAng-0.2),sampleRate,'hba','hba','h');


bodyMazeAng = [-xyz{t}(:,'bcom',[1,2]),...
               xyz{t}(:,'spine_upper',[1,2])-xyz{t}(:,'bcom',[1,2])];
bodyMazeAng = sq(bsxfun(@rdivide,bodyMazeAng,sqrt(sum(bodyMazeAng.^2,3))));
bodyMazeAng = cart2pol(bodyMazeAng(:,:,1),bodyMazeAng(:,:,2));
bodyMazeAng = circ_dist(bodyMazeAng(:,1),bodyMazeAng(:,2));
bodyMazeAng = MTADfet.encapsulate(Trials{tind},bodyMazeAng,sampleRate,'bma','bma','b');

hvec = xyz{t}(:,'nose',[1,2])-xyz{t}(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(rot(t)),-sin(rot(t));sin(rot(t)),cos(rot(t))],...
                 [2,3],...
                 [1,2]);
bvec = xyz{t}(:,'spine_upper',[1,2])-xyz{t}(:,'bcom',[1,2]);
bvec = sq(bsxfun(@rdivide,bvec,sqrt(sum(bvec.^2,3))));
bvec = cat(3,bvec,sq(bvec)*[0,-1;1,0]);
bvec = multiprod(bvec,...
                 [cos(rot(t)),-sin(rot(t));sin(rot(t)),cos(rot(t))],...
                 [2,3],...
                 [1,2]);

figure
t = 20;
u = 2;
ny = 5;
nx = 3;
[mxr,mxp] = pft{t}.maxRate(units{t}(u));        
pfstrj = MTADfet.encapsulate( ...
        Trials{t},...
        multiprod(bsxfun(@minus,...
                         mxp,...
                         sq(xyz{t}(:,'hcom',[1,2]))),...
                  hvec,2,[2,3]),...
        sampleRate,...
        'pfstrj','ptrj','t');


clf();
set(gcf(),'PaperType','A3');
set(gcf(),'PaperOrientation','landscape');
tspn = 1e4:2.6e4;
subplot(121);
    hold('on');
    scatter3(xyz{t}(tspn,'hcom',1),...
             xyz{t}(tspn,'hcom',2),300.*ones(size(tspn)),10,headBodyAng(tspn),'filled')
    plot(mxp(1),mxp(2),'or');
    plot(xyz{t}(tspn(1),'spine_middle',1),xyz{t}(tspn(1),'spine_middle',2),'^g');    
    plot(xyz{t}(tspn(end),'spine_middle',1),xyz{t}(tspn(end),'spine_middle',2),'^r');        
    plotSkeleton(Trials{t},xyz{t},tspn(1));        
    plotSkeleton(Trials{t},xyz{t},tspn(end));    
    xlim([-500,500]);
    ylim([-500,500]);    
    colormap('jet');    
    caxis([-1.2,1.2]);        
    title('Allocentric Coordinates')    
    daspect([1,1,1]);    
    grid('on');    
subplot(122);
    hold('on');
    scatter(pfstrj(tspn,2),pfstrj(tspn,1),10,headBodyAng(tspn),'filled')
    plot(pfstrj(tspn(1),2),pfstrj(tspn(1),1),'^g','MarkerSize',10,'MarkerFaceColor','g')
    plot(pfstrj(tspn(end),2),pfstrj(tspn(end),1),'^r','MarkerSize',10,'MarkerFaceColor','r')
    legend({'ego','start','stop'});
    xlim([-500,500]);    
    ylim([-500,500]);
    title('Egocentric Coordinates')
    colormap('jet');
    caxis([-1.2,1.2]);            
    daspect([1,1,1]);
    colorbar();
    grid('on');

    
    
figure();
ind = [Trials{t}.stc{'w+p+n&t',sampleRate}];
ind.cast('TimeSeries');
ind.resample(xyz{t});
% all
subplot2(ny,nx,1,1);
    indR = ind.data.*headBodyAng.data<-0.4;
    hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    colormap('jet');
    caxis([-1,1]);
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,1,2);
    indR = ind.data.*headBodyAng.data>-0.4&headBodyAng.data<0.4;
    hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    colormap('jet');
    caxis([-1,1]);
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,1,3);
    indL = ind.data.*headBodyAng.data>0.4;
    hist2([pfstrj(indL,2),pfstrj(indL,1)],linspace([-300,300,100]),linspace([-300,300,100]));    
    %scatter(pfstrj(indL,2),pfstrj(indL,1),1,headBodyAng(indL),'filled');
    colormap('jet');
    caxis([-1,1]);
    xlim([-300,300]);
    ylim([-300,300]);    
%hba < -0.4
subplot2(ny,nx,2,1); % 
    indR = ind.data & headBodyAng.data<-0.4 & bodyMazeAng.data < pi/4 & bodyMazeAng.data >-pi/4;
    hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('IN');
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,3,1); % OUT
    indL = ind.data & headBodyAng.data<-0.4 & (bodyMazeAng.data > pi*3/4 | bodyMazeAng.data < -pi*3/4);
    hist2([pfstrj(indL,2),pfstrj(indL,1)],linspace([-300,300,100]),linspace([-300,300,100]));    
    %scatter(pfstrj(indL,2),pfstrj(indL,1),1,headBodyAng(indL),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('OUT');    
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,4,1); % CW
    indR = ind.data & headBodyAng.data<-0.4 &  bodyMazeAng.data < -pi/4 & bodyMazeAng.data >-pi*3/4;
    hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('CW');
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,5,1); % CCW
    indL = ind.data & headBodyAng.data<-0.4 & (bodyMazeAng.data < pi*3/4 & bodyMazeAng.data > pi*1/4);
    hist2([pfstrj(indL,2),pfstrj(indL,1)],linspace([-300,300,100]),linspace([-300,300,100]));    
    %scatter(pfstrj(indL,2),pfstrj(indL,1),1,headBodyAng(indL),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('CCW');    
    xlim([-300,300]);
    ylim([-300,300]);    

%hba:0    
subplot2(ny,nx,2,2); % 
    indR = ind.data & headBodyAng.data<0.4&headBodyAng.data>-0.4 & bodyMazeAng.data < pi/4 & bodyMazeAng.data >-pi/4;
    hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('IN');
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,3,2); % OUT
    indL = ind.data & headBodyAng.data<0.4&headBodyAng.data>-0.4 & (bodyMazeAng.data > pi*3/4 | bodyMazeAng.data < -pi*3/4);
    hist2([pfstrj(indL,2),pfstrj(indL,1)],linspace([-300,300,100]),linspace([-300,300,100]));    
    %scatter(pfstrj(indL,2),pfstrj(indL,1),1,headBodyAng(indL),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('OUT');    
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,4,2); % CW
    indR = ind.data & headBodyAng.data<0.4&headBodyAng.data>-0.4 &  bodyMazeAng.data < -pi/4 & bodyMazeAng.data >-pi*3/4;
    hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('CW');
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,5,2); % CCW
    indL = ind.data & headBodyAng.data<0.4&headBodyAng.data>-0.4 & (bodyMazeAng.data < pi*3/4 & bodyMazeAng.data > pi*1/4);
    hist2([pfstrj(indL,2),pfstrj(indL,1)],linspace([-300,300,100]),linspace([-300,300,100]));    
    %scatter(pfstrj(indL,2),pfstrj(indL,1),1,headBodyAng(indL),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('CCW');    
    xlim([-300,300]);
    ylim([-300,300]);    
    
%hba > 0.4    
subplot2(ny,nx,2,3); % 
    indR = ind.data & headBodyAng.data>0.4& bodyMazeAng.data < pi/4 & bodyMazeAng.data >-pi/4;
    hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('IN');
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,3,3); % OUT
    indL = ind.data & headBodyAng.data>0.4& (bodyMazeAng.data > pi*3/4 | bodyMazeAng.data < -pi*3/4);
    hist2([pfstrj(indL,2),pfstrj(indL,1)],linspace([-300,300,100]),linspace([-300,300,100]));    
    %scatter(pfstrj(indL,2),pfstrj(indL,1),1,headBodyAng(indL),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('OUT');    
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,4,3); % CW
    indR = ind.data & headBodyAng.data>0.4&  bodyMazeAng.data < -pi/4 & bodyMazeAng.data >-pi*3/4;
    hist2([pfstrj(indR,2),pfstrj(indR,1)],linspace([-300,300,100]),linspace([-300,300,100]));
    %scatter(pfstrj(indR,2),pfstrj(indR,1),1,headBodyAng(indR),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('CW');
    xlim([-300,300]);
    ylim([-300,300]);    
subplot2(ny,nx,5,3); % CCW
    indL = ind.data & headBodyAng.data>0.4& (bodyMazeAng.data < pi*3/4 & bodyMazeAng.data > pi*1/4);
    hist2([pfstrj(indL,2),pfstrj(indL,1)],linspace([-300,300,100]),linspace([-300,300,100]));    
    %scatter(pfstrj(indL,2),pfstrj(indL,1),1,headBodyAng(indL),'filled');
    colormap('jet');
    caxis([-1,1]);
    title('CCW');    
    xlim([-300,300]);
    ylim([-300,300]);    



hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1.164, 44.873, 48.339, 4.868];
ny = 5;
nx = numel(phzBins)-1;
t = 18;
u = 2;


for u = 1:numel(units{t})
    clf();
    subplot2(ny,nx+1,1,1);
    plot(pft{t},units{t}(u),1,'text');
    title([Trials{t}.filebase,' - unit: ',num2str(units{t}(u))]);

    subplot2(ny,nx+1,2,1);    
    ind = ismember(spkvn.map,[t,units{t}(u)],'rows')               ... Unit
          & spkvn.stc(:,1)==1                                        ... Theta State
          & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9)  ... Bhv States
          & ismember(phzBinInd,[5,6]);                                          ... Theta phase
    scatter(spkvn.dst(ind),spkvn.bma(ind),10,chba(ind),'filled');
    colormap('jet');
    caxis([-1,1]);
% $$$ 
% $$$     subplot2(ny,nx+1,3,1);    
% $$$     [mxr,mxp] = pft{t}.maxRate(units{t}(u));        
% $$$     ind = [Trials{t}.stc{'w+p+n&t'}];
% $$$     scatter(xyz{t}(ind,'hcom',1)-mxp(1),...
% $$$             xyz{t}(ind,'hcom',2)-mxp(2),...
% $$$             1,...
% $$$             bodyMazeAng(ind),...
% $$$             'filled');
% $$$     xlim([-300,300]);
% $$$     ylim([-300,300]);
% $$$     Lines([],0,'k');
% $$$     Lines(0,[],'k');    

% $$$     subplot2(ny,nx+1,4,1);    
% $$$     pfstrj = MTADfet.encapsulate( ...
% $$$         Trials{t},...
% $$$         multiprod(bsxfun(@minus,...
% $$$                          mxp,...
% $$$                          sq(xyz{t}(:,'hcom',[1,2]))),...
% $$$                   hvec,2,[2,3]),...
% $$$         sampleRate,...
% $$$         'pfstrj','ptrj','t');
% $$$         
% $$$ 
% $$$     scatter(pfstrj(ind,1),...
% $$$             pfstrj(ind,2),...
% $$$             1,...
% $$$             bodyMazeAng(ind),...
% $$$             'filled');
% $$$     xlim([-300,300]);
% $$$     ylim([-300,300]);
% $$$     Lines([],0,'k');
% $$$     Lines(0,[],'k');    

    
    for x = 1:nx
% SELECT based on --> 
        ind = ismember(spkvn.map,[t,units{t}(u)],'rows')               ... Unit
              & spkvn.stc(:,1)==1                                        ... Theta State
              & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9)  ... Bhv States
              & phzBinInd == x;                                          ... Theta phase
% 
        subplot2(ny,nx+1,1,x+1);
            scatter(spkvn.ego(ind,2),   ...
                    spkvn.ego(ind,1)-25,      ...
                    10,                    ...
                    chba(ind),             ...        
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1.2,1.2]);
            colormap('jet');
            xlim([-300,300]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
% INBOUND 
        hpaInd = spkvn.bma < pi/4 & spkvn.bma >-pi/4;
        subplot2(ny,nx+1,2,x+1);
            scatter(spkvn.ego(ind&hpaInd,2),...
                    spkvn.ego(ind&hpaInd,1)-25,   ...
                    10,                    ...
                    chba(ind&hpaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&hpaInd&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1.2,1.2]);
            colormap('jet');
            xlim([-300,300]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
            if x==1,ylabel('IN');end
% OUTBOUND             
        hpaInd = spkvn.bma > pi*3/4 | spkvn.bma < -pi*3/4;            
        subplot2(ny,nx+1,3,x+1);
            scatter(spkvn.ego(ind&hpaInd,2),   ...
                    spkvn.ego(ind&hpaInd,1)-25,      ...
                    10,                    ...
                    chba(ind&hpaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&hpaInd&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1.2,1.2]);
            colormap('jet');
            xlim([-300,300]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
            if x==1,ylabel('OUT');end            
% CW 
        hpaInd = spkvn.bma > -pi*3/4 & spkvn.bma < -pi*1/4;                        
        subplot2(ny,nx+1,4,x+1);
            scatter(spkvn.ego(ind&hpaInd,2),   ...
                    spkvn.ego(ind&hpaInd,1)-25,      ...
                    10,                    ...
                    chba(ind&hpaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&hpaInd&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1.2,1.2]);
            colormap('jet');
            xlim([-300,300]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
% CCW 
        hpaInd = spkvn.bma < pi*3/4 & spkvn.bma > pi*1/4;            
        subplot2(ny,nx+1,5,x+1);
            scatter(spkvn.ego(ind&hpaInd,2),   ...
                    spkvn.ego(ind&hpaInd,1)-25,      ...
                    10,                    ...
                    chba(ind&hpaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&hpaInd&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1.2,1.2]);
            colormap('jet');
            xlim([-300,300]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
            if x==1,ylabel('CCW');end
            
    end

    waitforbuttonpress();
end



% HBA color
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1.164, 44.873, 48.339, 4.868];
ny = 5;
t = 18;
for u = 1:numel(units{t})
    clf();
    subplot2(ny,nx+1,1,1);
    plot(pft{t},units{t}(u),1,'text');
    title([Trials{t}.filebase,' - unit: ',num2str(units{t}(u))]);
    for x = 1:nx
% SELECT based on --> 
        ind = ismember(spkvn.map,[t,units{t}(u)],'rows')               ... Unit
              & spkvn.stc(:,1)==1                                        ... Theta State
              & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9)  ... Bhv States
              & phzBinInd == x;                                          ... Theta phase
% 
        subplot2(ny,nx+1,1,x+1);
            scatter(spkvn.ego(ind,1)-25,   ...
                    spkvn.ego(ind,2),      ...
                    10,                    ...
                    chba(ind),             ...        
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1,1]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
% E
        bmaInd = spkvn.bma < pi/4 & spkvn.bma >-pi/4;
        subplot2(ny,nx+1,2,x+1);
            scatter(spkvn.ego(ind&bmaInd,1)-25,   ...
                    spkvn.ego(ind&bmaInd,2),      ...
                    10,                    ...
                    chba(ind&bmaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&bmaInd&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1,1]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
        bmaInd = spkvn.bma > pi*3/4 | spkvn.bma < -pi*3/4;            
        subplot2(ny,nx+1,3,x+1);
            scatter(spkvn.ego(ind&bmaInd,1)-25,   ...
                    spkvn.ego(ind&bmaInd,2),      ...
                    10,                    ...
                    chba(ind&bmaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&bmaInd&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1,1]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
        bmaInd = spkvn.bma < pi*3/4 & spkvn.bma > pi*1/4;            
        subplot2(ny,nx+1,4,x+1);
            scatter(spkvn.ego(ind&bmaInd,1)-25,   ...
                    spkvn.ego(ind&bmaInd,2),      ...
                    10,                    ...
                    chba(ind&bmaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&bmaInd&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1,1]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
        bmaInd = spkvn.bma > -pi*3/4 & spkvn.bma < -pi*1/4;                        
        subplot2(ny,nx+1,5,x+1);
            scatter(spkvn.ego(ind&bmaInd,1)-25,   ...
                    spkvn.ego(ind&bmaInd,2),      ...
                    10,                    ...
                    chba(ind&bmaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.ego(ind&bmaInd&abs(spkvn.ego(:,1))<350,1)-25))));
            caxis([-1,1]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
    end
waitforbuttonpress();
end





hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1.164, 44.873, 48.339, 4.868];
ny = 3;
t = 20;
for u = 1:numel(units{t})
    clf();
    subplot2(ny,nx+1,1,1);
    plot(pft{t},units{t}(u),1,'text');
    title([Trials{t}.filebase,' - unit: ',num2str(units{t}(u))]);
    for x = 1:nx
% SELECT based on --> 
        ind = ismember(spkvn.map,[t,units{t}(u)],'rows')               ... Unit
              & spkvn.stc(:,1)==1                                        ... Theta State
              & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9)  ... Bhv States
              & phzBinInd == x;                                          ... Theta phase
% 
        subplot2(ny,nx+1,1,x+1);
            scatter(spkvn.bdy(ind,1)-25,   ...
                    spkvn.bdy(ind,2),      ...
                    10,                    ...
                    spkvn.hvl(ind),             ...        
                    'filled');
            title(num2str(round(median(spkvn.bdy(ind&abs(spkvn.bdy(:,1))<350,1)-25))));
            caxis([-10,10]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
% E
        bmaInd = spkvn.bma < pi/4 & spkvn.bma >-pi/4;
        subplot2(ny,nx+1,2,x+1);
            scatter(spkvn.bdy(ind&bmaInd,1)-25,   ...
                    spkvn.bdy(ind&bmaInd,2),      ...
                    10,                    ...
                    spkvn.hvl(ind&bmaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.bdy(ind&bmaInd&abs(spkvn.bdy(:,1))<350,1)-25))));
            caxis([-10,10]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
        bmaInd = spkvn.bma > pi*3/4 | spkvn.bma < -pi*3/4;            
        subplot2(ny,nx+1,3,x+1);
            scatter(spkvn.bdy(ind&bmaInd,1)-25,   ...
                    spkvn.bdy(ind&bmaInd,2),      ...
                    10,                    ...
                    spkvn.hvl(ind&bmaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.bdy(ind&bmaInd&abs(spkvn.bdy(:,1))<350,1)-25))));
            caxis([-10,10]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
        bmaInd = spkvn.bma < pi*3/4 & spkvn.bma > pi*1/4;            
        subplot2(ny,nx+1,4,x+1);
            scatter(spkvn.bdy(ind&bmaInd,1)-25,   ...
                    spkvn.bdy(ind&bmaInd,2),      ...
                    10,                    ...
                    spkvn.hvl(ind&bmaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.bdy(ind&bmaInd&abs(spkvn.bdy(:,1))<350,1)-25))));
            caxis([-10,10]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
        bmaInd = spkvn.bma > -pi*3/4 & spkvn.bma < -pi*1/4;                        
        subplot2(ny,nx+1,5,x+1);
            scatter(spkvn.bdy(ind&bmaInd,1)-25,   ...
                    spkvn.bdy(ind&bmaInd,2),      ...
                    10,                    ...
                    spkvn.hvl(ind&bmaInd),        ...                            
                    'filled');
            title(num2str(round(median(spkvn.bdy(ind&bmaInd&abs(spkvn.bdy(:,1))<350,1)-25))));
            caxis([-10,10]);
            colormap('jet');
            xlim([-300,400]);
            ylim([-300,300]);
            grid(gca(),'on');
            Lines([],0,'k');
            Lines(0,[],'k');
    end
waitforbuttonpress();
end





hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1.164, 44.873, 48.339, 9.868];
ny = numel(hvlBins)-1;
ny = 2;
nx = numel(phzBins)-1;
t = 20;
u = 3;
for u = 1:numel(units{t})
    clf();
    subplot2(ny,nx+1,1,1);
    plot(pft{t},units{t}(u),1,'text');
% $$$     [mxr,mxp] = pft{t}.maxRate(units{t}(u));        
% $$$     pfstrj = MTADfet.encapsulate( ...
% $$$         Trials{t},...
% $$$         multiprod(bsxfun(@minus,...
% $$$                          mxp,...
% $$$                          sq(xyz{t}(:,'hcom',[1,2]))),...
% $$$                   hvec,2,[2,3]),...
% $$$         sampleRate,...
% $$$         'pfstrj','ptrj','t');
% $$$ 

    
for x = 1:nx
    %for y = 1:ny
        %for v = 1:nv
        %subplot2(ny*nv,nx,y+3*(v-1),x);            
        subplot2(ny,nx+1,y,x+1);
% SELECT based on ->
        ind = ismember(spkvn.map,[t,units{t}(u)],'rows')               ...
              & spkvn.stc(:,1)==1                                       ... Theta State
              & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9) ... Bhv States
              & phzBinInd == x;                                          ... Theta phase
              %& hvfBinInd == y                                          ... forward Head speed
              %& hvlBinInd == y                                          ... lateral Head speed
              %& hbaBinInd == 2;                                           %  Head Body Angle
        
        scatter(spkvn.ego(ind,2),   ...
                spkvn.ego(ind,1)-25,      ...
                10,                    ...
                chba(ind),             ...        
                'filled');
        caxis([-1,1]);
        colormap('jet');
        xlim([-300,300]);
        ylim([-300,300]);
        grid(gca(),'on');
        Lines([],0,'k');
        Lines(0,[],'k');
        %end
    %end
end
waitforbuttonpress();
end

%ind = ismember(spkvn.map(:,1),[1,2,6,7])                     ... Session
ind = ismember(spkvn.map(:,1),[3,4,5,8:14])                     ... Session
      & spkvn.stc(:,1)==1                                       ... Theta State
      & phzBinInd == 8                                          ... Theta phase
      & hvlBinInd == 2                                          ... Theta phase
      & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9);  % Bhv States

figure,
scatter3(spkvn.ego(ind,1)-25,   ...
         spkvn.ego(ind,2)-8,   ...
         sin(chba(ind)),10,spkvn.hvl(ind),'filled');
xlim([-400,400]);
ylim([-400,400]);
zlim([-1,1]);
caxis([-60,60]);
colormap('jet');




figure,
hist2([spkvn.ego(ind,2)-8,sin(chba(ind))],linspace(-400,400,16),linspace(-1,1,16));






figure();
nx = numel(phzBins)-1;
for x = 1:nx
    ind = ismember(spkvn.map(:,1),[3,4,5,8:14])                     ... Session
          & spkvn.stc(:,1)==1                                       ... Theta State
          & phzBinInd == x                                          ... Theta phase
          & hvlBinInd == 2                                          ... Theta phase
          & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9);  % Bhv States
          
    ind(ind) = randn(sum(ind),1)>0;
    subplot(2,nx,x);
        hist2([spkvn.ego(ind,2)-8,chba(ind)],...
              linspace(-400,400,16),linspace(-pi/2,pi/2,16));
        [R,P] = corrcoef(spkvn.ego(ind,2)-8,chba(ind));
        title({['R: ', num2str(round(R(2),4))],[' P: ',num2str(P(2))]});    
    subplot(2,nx,x+nx);    
        hist2([spkvn.ego(ind,1)-25,chba(ind)],...
              linspace(-400,400,16),linspace(-pi/2,pi/2,16));
        [R,P] = corrcoef(spkvn.ego(ind,1)-8,chba(ind));
        title({['R: ', num2str(round(R(2),4))],[' P: ',num2str(P(2))]});
        
end



figure();
nx = numel(phzBins)-1;
for x = 1:nx
    ind = ismember(spkvn.map(:,1),[3,4,5,8:14])                     ... Session
          & spkvn.stc(:,1)==1                                       ... Theta State
          & phzBinInd == x                                          ... Theta phase
          & hvlBinInd == 2                                          ... Theta phase
          & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9);  % Bhv States
    ind(ind) = randn(sum(ind),1)>0;
    subplot(2,nx,x);
        hist2([spkvn.ego(ind,2)-8,chba(ind)],...
              linspace(-400,400,16),linspace(-1,1,16));
        [R,P] = corrcoef(spkvn.ego(ind,2)-8,chba(ind));
        title({['R: ', num2str(round(R(2),4))],[' P: ',num2str(P(2))]});    
    subplot(2,nx,x+nx);    
        hist2([spkvn.ego(ind,1)-25,sin(chba(ind))],...
              linspace(-400,400,16),linspace(-1,1,16));
        [R,P] = corrcoef(spkvn.ego(ind,1)-25,chba(ind));
        title({['R: ', num2str(round(R(2),4))],[' P: ',num2str(P(2))]});
end





[R,P] = corrcoef(atan2(spkvn.ego(ind,1)-25,spkvn.ego(ind,2)-8),sin(chba(ind)));
        

ny = numel(hvlBins)-1;
nx = numel(phzBins)-1;      
figure
for y = 1:ny;
for x = 1:nx;
ind = ismember(spkvn.map(:,1),[3,4,5,8:14])                     ... Session
      & spkvn.stc(:,1)==1                                       ... Theta State
      & phzBinInd == x                                          ... Theta phase
      & hvlBinInd == y                                          ... Theta phase
      & (spkvn.stc(:,7)==7|spkvn.stc(:,8)==8|spkvn.stc(:,9)==9);  % Bhv States
    subplot2(ny,nx,y,x);      
ce = [spkvn.ego(ind,2)-8,spkvn.ego(ind,1)-25];
ce = ce./sqrt(sum(ce.^2));
ce = atan2(ce(:,1),ce(:,2));
ce = ce+pi;
ce(ce>pi) = ce(ce>pi) - 2*pi;
%plot([chba(ind);chba(ind);chba(ind)],[ce-2*pi;ce;ce+2*pi],'.','MarkerSize',1);
hist2([[chba(ind);chba(ind);chba(ind)],[ce-2*pi;ce;ce+2*pi]],30,90)
Lines([],[-pi,pi],'k');
Lines(0,[],'k');
xlim([-pi,pi]);
end
end



% testing 
lims = {[-300,300],[-200,200],[-100,100]};
maskDimInd = cf(@(p,l) l(1) < p & p < l(2),  pft{t}.adata.bins([1:2,2]),lims);
mask = maskDimInd{1};
for dim = 2:numel(maskDimInd),
     mask = bsxfun(@and,mask,permute(maskDimInd{dim},[2:dim,1]));
end
figure,pcolor(pft{t}.adata.bins{1:2},double(mask(:,:,20))')