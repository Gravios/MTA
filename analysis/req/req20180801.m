%% phase precession during grooming and sitting
% spk
% phz
% drz

% find main region of sitting
ny = 2;
nx = 2;
figure();
i = 1;
subplot(ny,nx,i);i=i+1;
ind = stc{'m'}; % grooming periods
plot(xyz(ind,'nose',1),xyz(ind,'nose',2),'.');
xlim([-500,500]);
ylim([-500,500]);
subplot(ny,nx,i);i=i+1;
ind = stc{'s'}; % grooming periods
plot(xyz(ind,'nose',1),xyz(ind,'nose',2),'.');
xlim([-500,500]);
ylim([-500,500]);

% get units which have peak drz greater than 0.3 where y < -100
unitsGroom = unitSubset(mean(1-abs(drz(stc{'m'},:)))>0.2);
unitsSit   = unitSubset(mean(1-abs(drz(stc{'s'},:)))>0.2);



munits = unitSubset(1:69);
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,40,30];
ind = stc{'hpause+rear&t',xyz.sampleRate};
for u = 1:numel(munits),
    res = spk(munits(u));
    res = res(WithinRanges(res,ind.data));
    subplotfit(u,numel(munits));
    %hist2([repmat(ang(res,'pelvis_root','spine_upper',2),[2,1]),...           
% $$$     hist2([repmat(ang(res,'spine_middle','spine_upper',2),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi]],...
% $$$           linspace(-pi/2,pi/2,30),...
% $$$           linspace(-pi,pi*3,30));
% $$$     caxis([0,numel(res)/100]);
 
% $$$     plot(repmat(ang(res,'pelvis_root','nose',2),[2,1]),...
% $$$       ([phz(res,1);phz(res,1)+2*pi]),'.')
% $$$     %           flipud([phz(res,1);phz(res,1)+2*pi]),'.')
% $$$     xlim([-0.7,1.5]);

% $$$     %hist2([repmat(ang(res,'pelvis_root','spine_upper',2),[2,1]),...           
    hist2([repmat(ang(res,'spine_middle','hcom',2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi])],...
          linspace(-0.6,1.3,30),...
          linspace(-pi,pi*3,30));
    caxis([0,numel(res)/150]);
    

% $$$     plot(repmat(circ_dist(ang(res,'head_back','head_front',2),ang(res,'pelvis_root','spine_upper',2)),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi],'.')         
    %flipud([phz(res,2);phz(res,2)+2*pi]),'.')
% $$$     xlim([-1.5,1]);
    %hist2([repmat(ang(res,'pelvis_root','spine_upper',2),[2,1]),...           
% $$$     hist2([repmat(circ_dist(ang(res,'head_back','head_front',2),ang(res,'spine_middle','spine_upper',2)),[2,1]),...
% $$$            [phz(res,1);phz(res,1)+2*pi]],...
% $$$           linspace(-1.25,1.25,40),...
% $$$           linspace(-pi,pi*3,40));
% $$$     caxis([0,numel(res)/150]);
    
end

x = 60;
munits = unitSubset(x:x+9);
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,60,20];
% $$$ ind = stc{'rear&t',xyz.sampleRate};
% $$$ markers = {'spine_middle','spine_upper'}; closeupLims = [0.6,1.4];
ind = stc{'hpause&t',xyz.sampleRate};
markers = {'head_back','head_front'}; closeupLims = [-0.4,0.9];

for u = 1:numel(munits),
    res = spk(munits(u));
    res = res(WithinRanges(res,ind.data));
 
    % height 
% $$$     subplot2(3,10,1,u);
% $$$     plot(repmat(xyz(res,markers{2},3),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi]),'.');
% $$$     subplot2(3,10,2,u);
% $$$     plot(repmat(xyz(res,markers{2},3),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi]),'.');
% $$$     subplot2(3,10,3,u);    
% $$$     plot(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi],'.')
% $$$     
    % PITCH 
    subplot2(3,10,1,u);
    plot(repmat(ang(res,markers{1},markers{2},2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi]),'.');
    xlim(closeupLims);
    ylim([-pi,pi*3]);
 
    subplot2(3,10,2,u);
    plot(repmat(ang(res,markers{1},markers{2},2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi]),'.');
    xlim([-1,pi/2]);
    ylim([-pi,pi*3]);
    subplot2(3,10,3,u);    
    
    plot(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
           [phz(res,2);phz(res,2)+2*pi],'.')
    xlim([-1,1])
    ylim([-pi,pi*3]);
    
% $$$     % PITCH 
% $$$     subplot2(3,10,1,u);
% $$$     hist2([repmat(ang(res,'head_back','head_front',2),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi])],...
% $$$           linspace(-.4,0.6,20),...
% $$$           linspace(-pi,pi*3,20));
% $$$     caxis([0,numel(res)/80]);
% $$$ 
% $$$     subplot2(3,10,2,u);
% $$$     hist2([repmat(ang(res,'head_back','head_front',2),[2,1]),...
% $$$            ([phz(res,1);phz(res,1)+2*pi])],...
% $$$           linspace(-1,pi/2,20),...
% $$$           linspace(-pi,pi*3,20));
% $$$     caxis([0,numel(res)/50]);
% $$$     
% $$$     subplot2(3,10,3,u);    
% $$$     hist2([repmat(drz(res,munits(u)==unitSubset),[2,1]),...
% $$$            [phz(res,2);phz(res,2)+2*pi]],...
% $$$           linspace(-1,1,20),...
% $$$           linspace(-pi,pi*3,20));
% $$$     caxis([0,numel(res)/50]);
end

figure();
u = 9;
res = spk(munits(u));
res = res(WithinRanges(res,ind.data));
plot3(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
      repmat(ang(res,'head_back','head_front',2),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi]),'.');


x = 41;
munits = unitSubset(x:x+9);
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [1,1,60,10];
ind = stc{'hpause+rear&t',xyz.sampleRate};
for u = 1:numel(munits),
    res = spk(munits(u));
    res = res(WithinRanges(res,ind.data));
    
    % PITCH 
    subplot2(2,10,1,u);
    hist2([repmat(circ_dist(ang(res,'pelvis_root','spine_upper',2),ang(res,'head_back','head_front',2)),[2,1]),...
           ([phz(res,1);phz(res,1)+2*pi])],...
          linspace(-1.2,1,20),...
          linspace(-pi,pi*3,20));
    caxis([0,numel(res)/50]);
    
    subplot2(2,10,2,u);    
    hist2([repmat(drz(res,munits(u)==unitSubset),[2,1]),...
           [phz(res,2);phz(res,2)+2*pi]],...
          linspace(-1,1,20),...
          linspace(-pi,pi*3,20));
    caxis([0,numel(res)/50]);
end


munits = unitSubset(1:5);
figure();
ind = stc{'w',xyz.sampleRate};
for u = 1:numel(munits),
    res = spk(munits(u));
    res = res(WithinRanges(res,ind.data));
    %for c = 1:size(phz,2),
        subplot2(4,numel(unitsSit),c,u);
        hist2([repmat(drz(res,unitSubset==munits(u)),[2,1]),...
              [uiphz(res,c);uiphz(res,c)+2*pi]],...
              linspace(-1,1,30),...
              linspace(-pi,pi*3,30));
        caxis([0,numel(res)/100]);
    %end
end

    

figure(); ind = stc{'m'};
for u = 1:numel(unitsGroom),
    res = spk(unitsGroom(u));
    res = res(WithinRanges(res,ind.data));
    for c = 1:size(phz,2),
        subplot2(4,numel(unitsGroom),c,u);
        hist2([repmat(drz(res,unitSubset==unitsGroom(u)),[2,1]),...
              [phz(res,c);phz(res,c)+2*pi]],...
              linspace(-1,1,30),...
              linspace(-pi,pi*3,30));
        caxis([0,numel(res)/100]);        
    end
end


per = stc{'s-t'};

awper = per&ThreshCross(-lys.data,-1.7,2);
mergePeriods = ThreshCross((awper(2:end,1)-awper(1:end-1,2))<10,0.5,0);
while size(mergePeriods,1)~=0,
    awper.data(mergePeriods(1,1),:) = [awper(mergePeriods(1,1),1),awper(mergePeriods(1,2),2)];
    awper.data((mergePeriods(1,1)+1):mergePeriods(1,2),:) = [];
    mergePeriods = mergePeriods - [mergePeriods(1,2) - mergePeriods(1,1)];
    mergePeriods(1,:) = [];
end

        
        

figure();
ind = awper;
for u = 1:numel(unitsSit),
    res = spk(unitsSit(u));
    res = res(WithinRanges(res,ind.data));
    for c = 1:size(phz,2),
        subplot2(4,numel(unitsSit),c,u);
        rose(phz(res,c),30);
% $$$         hist2([repmat(drz(res,unitSubset==unitsSit(u)),[2,1]),...
% $$$               [phz(res,c);phz(res,c)+2*pi]],...
% $$$               linspace(-1,1,30),...
% $$$               linspace(-pi,pi*3,30));
% $$$         caxis([0,numel(res)/100]);
% $$$         plot(repmat(drz(res,unitSubset==unitsGroom(u)),[2,1]),...
% $$$               [phz(res,c);phz(res,c)+2*pi],...
% $$$              '.');
    end
end




sawper = awper.copy();
win = 150;
pdur = diff(awper.data,1,2);
sawper.data = awper(pdur>win,1);
sawper.data = cat(2,sawper.data,sawper.data+win);


tper = nan([size(phz,1),1]);
for p = 1:size(sawper,1),
    tper((sawper(p,1)):(sawper(p,1)+win-1)) = 1:win;
    %tper((sawper(p,1)-win):(sawper(p,1)+win)) = -win:win;
end


figure();
ind = sawper;
for u = 1:numel(munits),
    res = spk(munits(u));
    res = SelectPeriods(res,ind.data,'d',1,0);
    for c = 1:size(phz,2),
        subplot2(4,numel(munits),c,u);
% $$$         hist2([repmat(tper(res),[2,1]),...
% $$$               [phz(res,c);phz(res,c)+2*pi]],...
% $$$               linspace(0,win,30),...
% $$$               linspace(-pi,pi*3,30));
% $$$         caxis([0,numel(res)/170]);
        plot(repmat(tper(res),[2,1]),...
              [phz(res,c);phz(res,c)+2*pi],...
             '.');
    end
end

per = stc{'s-t'};
munits = unitsSit;

per = stc{'m'};
munits = unitsGroom;


figure();
ind = per; 
for u = 1:numel(munits),
    res = spk(munits(u));
    res = SelectPeriods(res,ind.data,'d',1,0);
    for c = 1:size(phz,2),
        subplot2(4,numel(munits),c,u);
        hist2([repmat(lrs(res),[2,1]),       ...   % Rad/LM theta power ratio
               [phz(res,c);phz(res,c)+2*pi]],...   % theta phase
                linspace(0.6,1.1,30),        ...   % Rad/LM theta power ratio bins
                linspace(-pi,pi*3,30));            % Theta phase bins
% $$$         hist2([repmat(lys(res),[2,1]),...
% $$$                [phz(res,c);phz(res,c)+2*pi]],...
% $$$               linspace(1.1,2.5,30),...              linspace(0.6,1.1,30),...
% $$$         linspace(-pi,pi*3,30));
        caxis([0,numel(res)/170]);
% $$$         plot(repmat(lys(res),[2,1]),...
% $$$               [phz(res,c);phz(res,c)+2*pi],...
% $$$              '.');
    end
end




% CCG between units 
T = [];
G = [];
ind = per; 
%ind = awper;
for u = 1:numel(munits),
    res = spk(munits(u));
    res = SelectPeriods(res,ind.data,'d',1,0);
    T = [T;res];
    G = [G;u*ones([numel(res),1])];
end
    
[mccg,tbins] = CCG(T,G,1,24,xyz.sampleRate,unique(G),'hz');
figure();
for u1 = 1:numel(munits),
    for u2 = 1:numel(munits),
       subplot2(numel(munits),numel(munits),u1,u2);
       plot(tbins,mccg(:,u1,u2));
       xlim(tbins([1,end]));
       ylim([0,15]);
    end
end



per = stc{'x'};
munits = unitsGroom;

figure();
ind = per; for u = 1:numel(munits),
    res = spk(munits(u));
    res = SelectPeriods(res,ind.data,'d',1,0);
    for c = 1:size(phz,2),
        subplot2(4,numel(munits),c,u);
        scatter(repmat(drz(res,munits(u)==unitSubset),[2,1]),...
              [phz(res,c);phz(res,c)+2*pi],...
                8,repmat(lrs(res),[2,1]),'filled');
    end
end


    

% units 44,59,73
sunits = [44,59,73];
figure()
for u = 1:numel(sunits)
    subplot2(1,numel(sunits),1,u);
    plot(pfs,sunits(u),'mean',true,[],true);
end

    

Trial = MTATrial.validate('jg05-20120309.cof.all');
lfp = Trial.load('lfp',80:89);

specArgsGamma = struct('nFFT',2^7,...
                       'Fs',  lfp.sampleRate,...
                       'WinLength',2^6,...
                       'nOverlap',2^6*0.875,...
                       'NW',3,...
                       'Detrend',[],...
                       'nTapers',[],...
                       'FreqRange',[50,250]);
   
[ysg,fsg,tsg] = fet_spec(Trial,lfp,[],[],[],specArgsGamma);




figure();
sp = [];
for s = 1:10,
    sp(s) = subplot(10,1,s);
    hold('on');
    imagesc(tsg,fsg,log10(ysg(:,:,s))');    
    plot((1:size(lfp,1))./lfp.sampleRate,nunity(lfp(:,s))*15+120);
    axis('xy');
    colormap(gca,'jet');

end
linkaxes(get(gcf,'Children'),'x');
ForAllSubplots('caxis([1,4.5])')

sp(5) = subplot(6,1,5);
hold('on');
plot([[1:size(vxy,1)]'./vxy.sampleRate,[1:size(vxy,1)]'./vxy.sampleRate],vxy.data)
for s = 1:numel(sunits)
    plot(./xyz.sampleRate,ones([numel(spk(sunits(s))),1])./s,'*m');
end
sp(6) = subplot(6,1,6);
haxSTS = plotSTC(stc,1,[],fliplr(states),fliplr('rbcgkmy'));
haxSTS.YTickLabel = fliplr(states);
xlabel('Time (s)');
linkaxes(get(gcf,'Children'),'x');



lysg = ysg.copy;
lysg.data = log10(ysg(:,:,4));


lysg.resample(xyz);

lysBins = linspace(1,2.4,30);
lysBinsInds = lysg.copy();
lysBinsInds.data = discretize(lys.data,lysBins);
lysBinsInds.data = lysBinsInds.data.*double(any(stcm==7,2));
lysBinsInds.data = lysBinsInds(spk(sunits(1)));


lvxy = vxy.copy();
lvxy.data(lvxy.data<1e-3) = 1e-3;

lvxy.data = log10(lvxy.data)
lysBins = linspace(-3,2,30);
lvxy.data = log10(lvxy.data)

lysBins = linspace(-3,2,30);
lysBinsInds = lvxy.copy();

lysBins = linspace(-1.2,1,25);
lysBinsInds = ang.copy();
lysBinsInds.data = circ_dist(ang(:,5,7,2),ang(:,3,4,2));

lysBinsInds.data = discretize(lysBinsInds(:,1),lysBins);
lysBinsInds.data = lysBinsInds.data.*double(any(stcm==2|stcm==3|stcm==4|stcm==5&stcm==8,2));
lysBinsInds.data = lysBinsInds(spk(u));
%p = -1;lysBinsInds.data = lysBinsInds(WithinRanges(phz(:,4),[p,p+0.2]));


nlysg = ysg.copy;
nlysg.data = log10(ysg(:,:,4));
nlysg.data(nlysg.data<prctile(nlysg.data(:),80)) = nan;
figure();
subplot(121);
lysXlysgSit = nan([numel(lysBins),size(lysg,2)]);
for i = 1:numel(lysBins),
    lysXlysgSit(i,:) = mean(nlysg(lysBinsInds.data==i,:),'omitnan');
end
imagesc(lysBins,fsg,lysXlysgSit'); 
axis xy;
colorbar();
colormap('jet');
subplot(122);
lysXlysgSit = nan([numel(lysBins),size(lysg,2)]);
for i = 1:numel(lysBins),
    lysXlysgSit(i,:) = std(nlysg(lysBinsInds.data==i,:),'omitnan');
end
imagesc(lysBins,fsg,lysXlysgSit'); 
axis xy;
colorbar();
colormap('jet');



unitSubset =
  Columns 1 through 13
     5    17    20    21    22    23    24    25    28    31    32    33    34
  Columns 14 through 26
    35    37    41    44    48    51    52    59    60    61    63    68    72
  Columns 27 through 39
    73    74    79    80    81    83    85    86    89    90    93    96    98
  Columns 40 through 52
   102   103   104   105   109   110   111   113   116   118   119   120   128
  Columns 53 through 65
   129   132   133   134   136   138   139   140   141   142   144   145   149
  Columns 66 through 69
   151   159   179   181
