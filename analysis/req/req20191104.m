


%%%<<< hello

xxx = 30

%%%>>>

phzOrder = [4,3,2,1,8,7,6,5];
l = 25;
sts = 1;
ind =   dpostI                      ...
        & duincI                      ...
        & dstcm(:,1)==1               ...
        & dstcm(:,1)~=2               ...        
        & logical(dstcm(:,stid(sts))) ...
        & ~any(logical(dstcm(:,[7,8])),2);

% JPDF( UINC, PUFR)
figure();
sp = tight_subplot(8,1,[0.01,0.01],[0.1],0.1);
for p = 1:8
    axes(sp(p));
    pind = ind & dphz == phzBinCenters(phzOrder(p));
    out = hist2([duinc(pind),...
                 dpufr(pind)],...
                 linspace(0,l,l+1),...
                 linspace(0,l,l+1));
    imagesc(linspace(0,l,l),...
            linspace(0,l,l),...
            log10(out)');
    axis('xy');
    
    line([0,l],[0,l]-2,'Color','m');
end
colormap jet


% JPDF( POST, UINC)
figure();
sp = tight_subplot(8,1,[0.01,0.01],[0.1],0.1);
for p = 1:8
    axes(sp(p));
    pind = ind & dphz == phzBinCenters(phzOrder(p));
    out = hist2([duinc(pind),...
                 log10(dpost(pind))],...
                 linspace(0,l,l+1),...
                 linspace(-3,0,l+1));
    imagesc(linspace(0,l,l),...
            linspace(-3,-0,l+1),...
            log10(out)');
    axis('xy');
    
    line([0,l],[0,0.02]-2,'Color','m');
end
colormap jet



mdiff = [];
sdiff = [];
ind =   dpostI                      ...
      & duincI                      ...
      & dstcm(:,1)==1               ...
      & dstcm(:,2)~=2               ...        
      & logical(dstcm(:,stid(sts))) ...
      & ~any(logical(dstcm(:,[7,8])),2);

ssmMask = false([size(dphz,1)./numel(phzBinCenters),1]);
ssmMask(1:64:end) = true;
ssmMask = reshape(repmat(ssmMask,[1,numel(phzBinCenters)])',[],1);


for p = 1:7,
    pind = find(ind & dphz == phzBinCenters(phzOrder(p)) & circshift(ssmMask,randsample(1:64,1)));    
    for o = p+1:8,
        oind = find(ind & dphz == phzBinCenters(phzOrder(o)) & circshift(ssmMask,randsample(1:64,1)));
        maxSamples = min([numel(oind),numel(pind)]);
        pinds = randsample(pind,maxSamples);
        oinds = randsample(oind,maxSamples);        
        
        mdiff(p,o) = mean(dpufr(pinds)-dpufr(oinds));
        sdiff(p,o) = std(dpufr(pinds)-dpufr(oinds));
    end
end

figure,
subplot(121);
imagesc(mdiff')
subplot(122);
imagesc(sdiff')



mpufr = [];
spufr = [];
for p = 1:8,    
    pind = find(ind & dphz == phzBinCenters(phzOrder(p)) & circshift(ssmMask,randsample(1:64,1)));
    for i = 1:100;
        mpufr(p,i) = mean(dpufr(randsample(pind,round(numel(pind)./4))));
        spufr(p,i) = std(dpufr(randsample(pind,round(numel(pind)./4))));
    end
end

muinc = [];
suinc = [];
mpost = [];
spost = [];
for p = 1:8,
    pind = find(ind & dphz == phzBinCenters(phzOrder(p)) & circshift(ssmMask,randsample(1:64,1)));
    for i = 1:100,
        mpost(p,i) = mean(dpost(randsample(pind,round(numel(pind)./4))));        
        spost(p,i) = std(dpost(randsample(pind,round(numel(pind)./4))));                
        muinc(p,i) = mean(duinc(randsample(pind,round(numel(pind)./4))));
        suinc(p,i) = std(duinc(randsample(pind,round(numel(pind)./4))));
    end
end

figure();
subplot(121);
hold('on');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mpufr,'k');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mpufr+spufr,'r');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mpufr-spufr,'r');
Lines(pi,[],'k');
Lines([],mean(mean(mpufr,2)),'r');
subplot(122);
hold('on');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,muinc,'k');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,muinc+suinc,'r');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,muinc-suinc,'r');
Lines(pi,[],'k');
Lines([],mean(mean(muinc,2)),'r');


figure();
subplot(121);
hold('on');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mpost,'k');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mpost+spost,'r');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mpost-spost,'r');
Lines(pi,[],'k');
Lines([],mean(mean(mpost,2)),'r');



figure();
subplot(121);
hold('on');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mean(mpufr,2));
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mean(mpufr,2)+mean(spufr,2),'r');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mean(mpufr,2)-mean(spufr,2),'r');
Lines(pi,[],'k');
Lines([],mean(mean(mpufr,2)),'r');
subplot(122);
hold('on');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mean(muinc,2));
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mean(muinc,2)+mean(suinc,2),'r');
plot(phzBinCenters(phzOrder)+double(phzBinCenters(phzOrder)<0).*2*pi,mean(muinc,2)-mean(suinc,2),'r');
Lines(pi,[],'k');
Lines([],mean(mean(muinc,2)),'r');








hbangBinEdges = linspace(-pi/2,pi/2,11);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
hbangBinInd = discretize(dhbang+0.4-0.8*double(~ismember(dtind,[3,4,5])),hbangBinEdges);

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI;        

figure();
for f = 1:2,    
    for b = 1:10,
        indb = ind & hbangBinInd==b;
        subplot2(2,10,f,b);
        hold('on');
        hist2([ferr{f}(indb),                   ...
              dphz(indb)],                      ...
              ferrorBinEdges{f},              ...
              phzBins);
        hist2([ferr{f}(indb),                   ...
              dphz(indb)+2*pi],                      ...
              ferrorBinEdges{f},              ...
              phzBins+2*pi);
        axis('tight');
        title(num2str(hbangBinCenters(b)));
    end
end




hvangBinEdges = linspace(-0.6,0.6,11);
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
hvangBinInd = discretize(dhvang,hvangBinEdges);

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI;% ...
%        & ismember(dtind,[3,4,5]);        

eds = 6;
hbangBinEdges = linspace(-pi/2,pi/2,eds);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
hbangBinInd = discretize(dhbang+0.4-0.8*double(~ismember(dtind,[3,4,5])),hbangBinEdges);

xyvBinEdges = linspace(-0.4,1.8,eds);
xyvBinCenters = mean([xyvBinEdges(2:end); xyvBinEdges(1:end-1)]);
xyvBinInd = discretize(log10(dxyvel),xyvBinEdges);

hvangBinEdges = linspace(-0.15,0.15,eds);
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
hvangBinInd = discretize(dhvang,hvangBinEdges);
hvangFSBinInd = discretize(circshift(dhvang,-256),hvangBinEdges);
hvangPSBinInd = discretize(circshift(dhvang,256),hvangBinEdges);

hvangBinEdges = linspace(-0.15,0.15,eds);
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
hvangBinInd = discretize(dhvang,hvangBinEdges);


g = 1;
figure();
for f = 1:eds-1,    
    for b = 1:eds-1,
        indb = ind & hvangFSBinInd==b & hbangBinInd==f;
        subplot2(eds-1,eds-1,f,b);
        hold('on');
        out = hist2([ferr{g}(indb),                   ...
              dphz(indb)],                      ...
              ferrorBinEdges{g},              ...
              phzBins);
        imagesc(ferrorBinEdges{g},              ...
               phzBins,...
               imgaussfilt(out,[1,0.5])');
        imagesc(ferrorBinEdges{g},              ...
               phzBins+2*pi,...
               imgaussfilt(out,[1,0.5])');
        axis('tight');
        axis('xy');
        title(['av:',num2str(hvangBinCenters(b)),' hba:',num2str(hbangBinCenters(f))]);
    end
end
ForAllSubplots('xlim([-250,250]);');
ForAllSubplots('colormap(''jet'');');
ForAllSubplots('Lines([],pi/2,''k'');');
ForAllSubplots('Lines(0,[],''k'');');


figure();
for b = 1:4
ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[b+2])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI;        
subplot(2,2,b);
hist2([dhvang(ind),dhbang(ind)],linspace(-0.2,0.2,50),linspace(-pi/2,pi/2,50));
end




eds = 10;
hrvfBinEdges = linspace(-30,80,eds);
hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
hrvfBinInd = discretize(-dhrvf,hrvfBinEdges);

hrvlBinEdges = linspace(-60,60,eds);
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
hrvlBinInd = discretize(dhrvl,hrvlBinEdges);

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI;        


g = 2;
figure();
sax = reshape(tight_subplot(eds-1,eds-1,0.001,0,0)',[1,1].*eds-1)';

for f = 1:eds-1,    
    for b = 1:eds-1,
        indb = ind & hrvfBinInd==f & hrvlBinInd==b;
        axes(sax(f,b));
        %subplot2(eds-1,eds-1,f,b);
        hold('on');
        out = hist2([ferr{g}(indb),                   ...
              dphz(indb)],                      ...
              ferrorBinEdges{g},              ...
              phzBins);
        imagesc(ferrorBinEdges{g},              ...
               phzBins,...
               out');
        imagesc(ferrorBinEdges{g},              ...
               phzBins+2*pi,...
               out');
        axis('tight');
        axis('xy');
        hold(sax(f,b),'on');
        text(sax(f,b),...
             -220,...
             -pi,...
             ['lat:',num2str(round(hrvlBinCenters(b),2)), ...
              ' fwd:',num2str(round(hrvfBinCenters(f),2))],...
             'Color','w',...
             'Rotation',90);
    end
end
ForAllSubplots('xlim([-250,250])');
ForAllSubplots('colormap(''jet'')');
ForAllSubplots('Lines([],pi/2,''k'')');
ForAllSubplots('Lines(0,[],''k'')');


eds = 6;
hrvfBinEdges = linspace(-50,75,eds);
hrvfBinEdges = [-40,-5, 5, 15, 25, 75];
hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
hrvfBinInd = discretize(-dhrvf,hrvfBinEdges);

hvangBinEdges = linspace(-0.15,0.15,eds);
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
hvangBinInd = discretize(dhvang,hvangBinEdges);


ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[4])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI;

g = 1;
figure();
for f = 1:eds-1,    
    for b = 1:eds-1,
        indb = ind & hrvfBinInd==f & hvangBinInd==b;
        subplot2(eds-1,eds-1,f,b);
        hold('on');
        out = hist2([ferr{g}(indb),                   ...
              dphz(indb)],                      ...
              ferrorBinEdges{g},              ...
              phzBins);
        imagesc(ferrorBinEdges{g},              ...
               phzBins,...
               out');
        imagesc(ferrorBinEdges{g},              ...
               phzBins+2*pi,...
               out');
        axis('tight');
        axis('xy');
        title(['avh:',num2str(hvangBinCenters(b)),' fwd:',num2str(hrvfBinCenters(f))]);
    end
end
ForAllSubplots('xlim([-250,250])');
ForAllSubplots('colormap(''jet'')');
ForAllSubplots('Lines([],pi/2,''k'')');
ForAllSubplots('Lines(0,[],''k'')');







figure();hist(dhvang(ind),hvangBinEdges);
figure();hist(dhbang(ind)+0.4-0.8*double(~ismember(dtind(ind),[3,4,5])),linspace(-pi*2/3,pi*2/3,400));

%hrvlBinEdges = linspace(-60,60,eds);
hrvlBinEdges = [-60,-20,-15,-8,-3,3,8,15,20,60];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edx = numel(hrvlBinCenters);
%hrvlBinInd = discretize(circshift(-dhrvl,-25*8),hrvlBinEdges);
hrvlBinInd = discretize(circshift(-dhrvl,0),hrvlBinEdges);
%hrvlBinInd = discretize(circshift(-dhrvl,25*8),hrvlBinEdges);

hrvfBinEdges = linspace(-60,60,eds);
hrvfBinEdges = [-30,-5,0,5,10,20,30,50,80];
hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
edy = numel(hrvfBinCenters);
hrvfBinInd = discretize(-dhrvf,hrvfBinEdges);



%hvangBinEdges = linspace(-0.15,0.15,eds);
%hvangBinEdges = [-0.3,-0.06,-0.04,-0.024,-0.012,0.012,0.024,0.04,0.06,0.3];
hvangBinEdges = [ 0, 0.009, 0.018, 0.03, 0.05, 0.17];
%hvangBinEdges = linspace(-0.1,0.1,8);
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edy = numel(hvangBinCenters);
%hvangBinInd = discretize(circshift(dhvang,-25*8),hvangBinEdges);
hvangBinInd = discretize(abs(dhvang),hvangBinEdges);
%hvangBinInd = discretize(circshift(dhvang,25*8),hvangBinEdges);

hbangBinEdges = linspace(-pi/2,pi/2,11);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.25-0.5*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);




ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<360; % ...
%        & -dhrvf>=-10;




nout = {};
%for s = 1:5,
%shift = 0;
shift = 0;
e = 2;
figure();
sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
s = 1;
tind = ind;% & ismember(hrollBinInd,5:7);
set(gcf,'Name','all');
% $$$ tind = ind & ismember(hrollBinInd,1:3);
% $$$ set(gcf,'Name','Roll Left');
tfet = ferr{e}(tind);
tphz = dphz(tind);
txi = hbangBinInd(tind);
tyi = hvangBinInd(tind);
%tyi = hrollBinInd(tind);
for f = 1:edy,    
    for b = 1:edx
        axes(sax(f,b));       
        
        %indb = ind & hvangBinInd==f & hbangBinInd==b & smMask;%  & ~(dhrvl> -5 &dhrvl<5);
        %indb = ind & hrvlBinInd==b & hvangBinInd==f & ismember(hbangBinInd,9);
        %indb = ind & hvangBinInd==f & hbangBinInd==b & ismember(hrollBinInd,1:3); % left
        %indb = ind & hvangBinInd==f & hbangBinInd==b & ismember(hrollBinInd,5:7); %right
        indb = tyi==f & txi==b; %right        
        %indb = ind & hvangBinInd==f & hbangBinInd==b;
        %indb = ind & hrvlBinInd==f & hvangBinInd==b;        
        %subplot2(eds-1,eds-1,f,b);
        hold('on');
        out = hist2([tfet(indb),           ...
                     tphz(indb)],                    ...
                    ferrorBinEdges{e},              ...
                    phzBins);

        imagesc(ferrorBinEdges{e},            ...
                phzBins,                      ...
                imgaussfilt(out,[2,0.1])');                

        imagesc(ferrorBinEdges{e},            ...
                phzBins+2*pi,                 ...
                imgaussfilt(out,[2,0.1])');

        axis('tight');
        axis('xy');
% $$$         text(sax(f,b),...
% $$$              -220,...
% $$$              -pi,...
% $$$              ['hba:',num2str(round(hbangBinCenters(b),2)),...
% $$$               ' hva:',num2str(round(hrollBinCenters(f),2))],...
% $$$              'Color','w',...
% $$$              'Rotation',90);
        text(sax(f,b),...
             -220,...
             -pi,...
             ['hba:',num2str(round(hbangBinCenters(b),2)),...
              ' hva:',num2str(round(hvangBinCenters(f),2))],...
             'Color','w',...
             'Rotation',90);
        %outn = out./sum(out(:));
        %nout{s}(f,b) = sum(sum(outn.*log2(outn./(sum(outn,2,'omitnan')*sum(outn,'omitnan'))),'omitnan'),'omitnan');
    end
end
ForAllSubplots('xlim([-200,200]);');
ForAllSubplots('colormap(''jet'');');
ForAllSubplots('Lines([],pi,''k'');');
ForAllSubplots('Lines(0,[],''k'');');
%ForAllSubplots('caxis(caxis()./1.5);');
%set(gcf,'Name','all');








hvangBinEdges = [-0.3,-0.06,-0.04,-0.024,-0.012,0.012,0.024,0.04,0.06,0.3];
%hvangBinEdges = [-0.06,-0.04,-0.024,-0.018,-0.012,-0.006,0.006,0.012,0.018,0.024,0.04,0.06];
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edy = numel(hvangBinCenters);
hvangBinInd = discretize(circshift(dhvang,0),hvangBinEdges);

hrvlBinEdges = [-60,-20,-15,-8,-3,3,8,15,20,60];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edx = numel(hrvlBinCenters);
%edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(circshift(-dhrvl,0),hrvlBinEdges);


hbangBinEdges = linspace(-pi/2,pi/2,13);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,5])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380; % ...

shift = 0;
e = 2;
figure();
sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
s = 1;
for f = 1:edy,    for b = 1:edx
        axes(sax(f,b));       
        indb = ind & hrvlBinInd==b & hvangBinInd==f & ismember(hbangBinInd,8:10);
        %indb = ind & hvangBinInd==f & hbangBinInd==b;
        %indb = ind & hrvlBinInd==f & hbangBinInd==b;
        hold('on');
        out = hist2([ferr{e}(indb),dphz(indb)],  ferrorBinEdges{e},  phzBins);
        imagesc(ferrorBinEdges{e},    phzBins,      imgaussfilt(out,[2,0.1])');                
        imagesc(ferrorBinEdges{e},    phzBins+2*pi, imgaussfilt(out,[2,0.1])');
        axis('tight');        axis('xy');
% $$$         text(sax(f,b), -220, -pi,...
% $$$              ['hba:',num2str(round(hbangBinCenters(b),2)),...
% $$$               ' hva:',num2str(round(hvangBinCenters(f),2))],...
% $$$              'Color','w',...
% $$$              'Rotation',90);
        text(sax(f,b), -220, -pi,...
             ['hrl:',num2str(round(hrvlBinCenters(b),2)),...
              ' hva:',num2str(round(hvangBinCenters(f),2))],...
             'Color','w',...
             'Rotation',90);
% $$$         text(sax(f,b), -220, -pi,...
% $$$              ['hba:',num2str(round(hbangBinCenters(b),2)),...
% $$$               ' hrl:',num2str(round(hrvlBinCenters(f),2))],...
% $$$              'Color','w',...
% $$$              'Rotation',90);
end;    end
ForAllSubplots('xlim([-250,250]);');    ForAllSubplots('colormap(''jet'');');
ForAllSubplots('Lines([],0,''k'');');  ForAllSubplots('Lines([],pi.*2,''k'');');  
ForAllSubplots('Lines([],pi,''r'');');  ForAllSubplots('Lines(0,[],''k'');');








hvangBinEdges = [-0.3,-0.06,-0.04,-0.024,-0.012,0.012,0.024,0.04,0.06,0.3];
hvangBinEdges = [-0.06,-0.04,-0.024,-0.018,-0.012,-0.006,0.006,0.012,0.018,0.024,0.04,0.06];
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edy = numel(hvangBinCenters);
edx = numel(hvangBinCenters);
hvangBinInd = discretize(circshift(dhvang,0),hvangBinEdges);

hbangBinEdges = linspace(-pi/2,pi/2,12);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.25-0.5*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);

hrollBinEdges = linspace(-0.5,0.5,8);
hrollBinCenters = mean([hrollBinEdges(2:end); hrollBinEdges(1:end-1)]);
%edx = numel(hrollBinCenters);
edy = numel(hrollBinCenters);
hrollBinInd = discretize((dhroll-0.18-0.2*double(~ismember(dtind,[3,4,5]))),hrollBinEdges);
cdhroll = (dhroll-0.18-0.2*double(~ismember(dtind,[3,4,5])));


brefBinEdges = linspace(-5,5,11);
brefBinCenters = mean([brefBinEdges(2:end); brefBinEdges(1:end-1)]);
edx = numel(brefBinCenters);
brefBinInd = discretize(dbreff,brefBinEdges);

brefdBinEdges = linspace(-0.6,0.6,11);
brefdBinCenters = mean([brefdBinEdges(2:end); brefdBinEdges(1:end-1)]);
edy = numel(brefdBinCenters);
brefdBinInd = discretize(dbrefd,brefdBinEdges);

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380; % ...





shift = 0;
e = 2;
figure();
sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
s = 1;
for f = 1:edy,    for b = 1:edx
        axes(sax(f,b));       
        %indb = ind & hrvlBinInd==b & hvangBinInd==f & ismember(hbangBinInd,9);
        %indb = ind & hvangBinInd==f & hbangBinInd==b;
        indb = ind & hrollBinInd==f & hbangBinInd==b;
        %indb = ind & hrollBinInd==f & hbangBinInd==b;% & ismember(hvangBinInd,[1:4,8:11]);        
        %indb = ind & hrollBinInd==f & hvangBinInd==b & ismember(hbangBinInd,8:11);        
        %indb = ind & hrollBinInd==f & hvangBinInd==b & ismember(hbangBinInd,1:6);                
        %indb = ind & hrollBinInd==b & hvangBinInd==f & ismember(hbangBinInd,4:6);                
        hold('on');
        out = hist2([ferr{e}(indb),dphz(indb)],  ferrorBinEdges{e},  phzBins);
        imagesc(ferrorBinEdges{e},    phzBins,      imgaussfilt(out,[2,0.1])');                
        imagesc(ferrorBinEdges{e},    phzBins+2*pi, imgaussfilt(out,[2,0.1])');
        axis('tight');        axis('xy');
% $$$         text(sax(f,b), -220, -pi,...
% $$$              ['hro:',num2str(round(hrollBinCenters(f),2)),...
% $$$               ' hva:',num2str(round(hvangBinCenters(b),2))],...
% $$$              'Color','w',...
% $$$              'Rotation',90);
        text(sax(f,b), -220, -pi,...
             ['hba:',num2str(round(hbangBinCenters(b),2)),...
              ' hro:',num2str(round(hrollBinCenters(f),2))],...
             'Color','w',...
             'Rotation',90);
end;    end
ForAllSubplots('ylim([0,2*pi]);'); 
ForAllSubplots('xlim([-250,250]);'); 
ForAllSubplots('colormap(''jet'');');
ForAllSubplots('Lines(0,[],''k'');');
ForAllSubplots('Lines([],0,''k'');');  
ForAllSubplots('Lines([],pi,''m'');');
ForAllSubplots('Lines([],2*pi,''k'');');




figure();
sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
for f = 1:edy,    for b = 1:edx
        axes(sax(f,b));
        indb = ind & brefBinInd==b & brefdBinInd==f & ismember(hbangBinInd,7:9);
        hold('on');
        out = hist2([ferr{e}(indb),dphz(indb)],  ferrorBinEdges{e},  phzBins);
        imagesc(ferrorBinEdges{e},    phzBins,      imgaussfilt(out,[2,0.1])');                
        imagesc(ferrorBinEdges{e},    phzBins+2*pi, imgaussfilt(out,[2,0.1])');
        axis('tight');        axis('xy');
end;    end;
ForAllSubplots('xlim([-250,250]);');    ForAllSubplots('colormap(''jet'');');
ForAllSubplots('Lines([],pi,''k'');');  ForAllSubplots('Lines(0,[],''k'');');









ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380; % ...


figure
for p = 1:8,
    for b = 1:11,
        subplot2(8,11,p,b);
        pind = mod(p+4,8)+double((p+4)==8);
        indb = ind & dphz==phzBinCenters(pind)& hbangBinInd==b;
        out = hist2([ferr{2}(indb),           ...
                     ferr{1}(indb)],                    ...
                    ferrorBinEdges{2},              ...
                    ferrorBinEdges{1});
        
        %out = out./sum(out(:));
        imagesc(ferrorBinCenters{2},            ...
                ferrorBinCenters{1},   ...
                imgaussfilt(out,[2,2])');      
        title(num2str(circ_rad2ang(phzBinCenters(pind)+2*pi.*double(phzBinCenters(pind)<0))));
        xlim([-350,350]);
        ylim([-350,350]);
        Lines([],0,'w');
        Lines(0,[],'w');
        colormap(gca(),'jet');
        %caxis([0,500]);
        axis('xy');
    end
end





b=9;f = 9;
indb = ind & hrvlBinInd==f & hvangBinInd==b & (dhbang < -0.5 & dhbang > -0.9) & smMask;
figure,
subplot(121);
plot([ferr{1}(indb&dphz==phzBinCenters(3)),ferr{1}(indb&dphz==phzBinCenters(7))],...
     [ferr{2}(indb&dphz==phzBinCenters(3)),ferr{2}(indb&dphz==phzBinCenters(7))],'.')
subplot(122);
plot([ferr{1}(indb&dphz==phzBinCenters(3)),ferr{1}(indb&dphz==phzBinCenters(7))]',...
     [ferr{2}(indb&dphz==phzBinCenters(3)),ferr{2}(indb&dphz==phzBinCenters(7))]')


 


% ABS angvel and hbang

eds = 6;
ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380; % ...
%        & -dhrvf>=-10;


%hvangBinEdges = linspace(-0.15,0.15,eds);

hvangBinEdges =prctile(abs(dhvang(ind&abs(dhvang)<0.3)),linspace(0.001,99.999,eds));


% ABSOLUTE value angles
hvangBinEdges = [ 0, 0.012, 0.024,0.04,0.060 0.3];
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edy = numel(hvangBinCenters);
hvangBinInd = discretize(abs(dhvang),hvangBinEdges);

eds = 9;
hbangBinEdges = linspace(0,pi/2,eds);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(abs(dhbang+0.4-0.8*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);



uquad = (sign(dhbang).*sign(dhvang))<0;
usign = sum([sign(dhbang)==1&sign(dhvang)==-1,sign(dhbang)==-1&sign(dhvang)==1]*[1;-1],2);


uquad = (sign(dhbang).*sign(dhvang))>0;
usign = sum(double([sign(dhbang)==1&sign(dhvang)==1,sign(dhbang)==-1&sign(dhvang)==-1])*[-1;1],2);

dirPath = '/storage/share/Projects/BehaviorPlaceCode/decode/theta_resolved';

%sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';

stsLbls = {    'all','high','low','loc','pause'};
stsGrps = {[3,4,5,6], [3,4],[5,6],[3,5],[4,6]};

shifts = 0:8:2^8;

for s = 1:numel(stsLbls),
    ind =   logical(dstcm(:,1))                             ...
            & any(logical(dstcm(:,stsGrps{s})),2)           ...
            & dpostI                                        ...
            & duincI ...
            & dhdist<380; % ...
    
    for q = 1:6
        switch q
          case 1
            uquad = sign(dhbang)==-1&sign(dhvang)==1;
            usign = ones(size(dhbang));
            uLblq = 'LRQ';
            uDesc = {'Head     Direction, Motion:    ',...
                     '            CW       CCW       '};
          case 2        
            uquad = sign(dhbang)==1&sign(dhvang)==-1;
            usign = -ones(size(dhbang));
            uLblq = 'ULQ';
            uDesc = {'Head     Direction, Motion:    ',...
                     '            CCW       CW       '};
          case 3
            uquad = sign(dhbang)==-1&sign(dhvang)==-1;
            usign = -ones(size(dhbang));
            uLblq = 'LLQ';
            uDesc = {'Head     Direction, Motion:    ',...
                     '             CW       CW       '};
          case 4
            uquad = sign(dhbang)==1&sign(dhvang)==1;
            usign = ones(size(dhbang));
            uLblq = 'URQ';
            uDesc = {'Head     Direction, Motion:    ',...
                     '             CCW       CCW       '};
          case 5
            uquad = (sign(dhbang).*sign(dhvang))>0;
            usign = sum(double([sign(dhbang)==1&sign(dhvang)==1,sign(dhbang)==-1&sign(dhvang)==-1])*[-1;1],2);
            uLblq = 'ALL_SD'
            uDesc = {'Head     Direction, Motion:    ',...
                     '            CCW      CCW       ',...
                     'flip lat    CW       CW        '};
          case 6
            uquad = (sign(dhbang).*sign(dhvang))<0;
            usign = sum(double([sign(dhbang)==1&sign(dhvang)==1,sign(dhbang)==-1&sign(dhvang)==-1])*[-1;1],2);            
            %usign = sum([sign(dhbang)==1&sign(dhvang)==-1,sign(dhbang)==-1&sign(dhvang)==1]*[1;-1],2);
            uLblq = 'ALL_DD'
            uDesc = {'Head     Direction, Motion:    ',...
                     '            CW       CCW      ',...
                     'flip lat    CCW      CW       '};
        end
        
        uLblq = [uLblq, '_', stsLbls{s}];
        uDesc = [uDesc,{['States: ',stsLbls{s}]}];
        
        for g = fliplr(1:2),
            switch g
              case 1
                fileBase = ['egocentricHead_phasePrecession_decodedLon', '_', uLblq];
                usign = ones(size(dhbang));
              case 2
                fileBase = ['egocentricHead_phasePrecession_decodedLat', '_', uLblq];
            end
            

            [hfig, fig, fax, sax] = set_figure_layout(figure(666061),'A4','landscape',[],3,3,0.1,0.1);        
            for f = 1:edy,    
                for b = 1:edx
                    sax(f,b) = axes('Units','centimeters',                                         ...
                                    'Position',[fig.page.xpos(b)+xOffSet,                          ...
                                        fig.page.ypos(edy-f+1)+yOffSet,                          ...
                                        fig.subplot.width,     ...
                                        fig.subplot.height],   ...
                                    'FontSize', 8,                                                 ...
                                    'LineWidth',1);
                    
                    
                    indb = ind & hvangBinInd==f & hbangBinInd==b & uquad;
                    %indb = ind & hvangBinInd==b & hbangBinInd==f & uquad & circshift(smMask, randsample(shifts,1));
                    indb(indb) = randn([sum(indb),1]) > 0;
                    hold('on');

                    out = hist2([ferr{g}(indb).*usign(indb),        ... 
                                 dphz(indb)],                       ...
                                ferrorBinEdges{g},                  ...
                                phzBins);
                    %out = bsxfun(@rdivide,out,sum(out));
                    imagesc(ferrorBinEdges{g},                      ...
                            phzBins,                                ...
                            imgaussfilt(out,[2,0.1])');                
                    imagesc(ferrorBinEdges{g},                      ...
                            phzBins+2*pi,                           ...
                            imgaussfilt(out,[2,0.1])');
                    %imgaussfilt(out,[2,0.1])');
                    axis('tight');
                    axis('xy');
                    text(sax(f,b),...
                         -220,...
                         -pi,...
                         ['avh:',num2str(round(hvangBinCenters(f),2)), ...
                          ' hba:',num2str(round(hbangBinCenters(b),2))],...
                         'Color','w',...
                         'Rotation',90);
                    Lines([],pi/2,'k');
                    Lines( 0,  [],'k');
                    sax(f,b).YTickLabels = {};
                    sax(f,b).XTickLabels = {};
                end
            end
            
            axes(fax);
            % PRINT Discription
            text(fig.page.marginLeft,sum(fig.page.ypos(1)+fig.subplot.height)+1,[{fileBase},uDesc]);

            af(@(h) xlim(h,[-300,300]), sax);
            af(@(h) colormap(h,'jet'),  sax);

            % SAVE Figure as {png, pdf, eps}
            filePath = fullfile(dirPath,fileBase);
            print(hfig,'-dpng', [filePath,'.png']);
            print(hfig,'-depsc',[filePath,'.eps']);
        end    
    end
end




% ABSOLUTE value angles
hvangBinEdges = [ 0, 0.012, 0.024,0.04,0.060 0.3];
hvangBinEdges = [ 0, 0.012, 0.024,0.04,0.060 0.3];
hvangBinEdges = [ 0, 0.009, 0.018, 0.03, 0.05, 0.17];
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edy = numel(hvangBinCenters);
hvangBinInd = discretize(abs(dhvang),hvangBinEdges);

eds = 7;
hbangBinEdges = linspace(0,pi/2,eds);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(abs(dhbang+0.25-0.5*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);


numIter = 100;
pxavn = zeros( [numel( errorBinCenters ), ...
               numel( phzBinCenters   ), ...
               edx,                      ...
               edy,                      ...
               numel(stsGrps),           ...
               2,                        ...
               2,                        ...
               numIter,                  ...               
             ]                           ...
           );
% pxav( x, y, hba, hav, sts, fet, quad, itr)
shifts = 0:8:2^8;
for sts = 1:numel(stsGrps),
    ind =   logical(dstcm(:,1))                             ...
            & any(dstcm(:,stsGrps{sts}),2)                  ...
            & dpostI                                        ...
            & duincI ...
            & dhdist<380;
    disp(['[Info] State: ',stsLbls{sts}]);
    for g = 1:2,        
        for e = fliplr(1:2),
            switch g
              case 1
                uquad = (sign(dhbang).*sign(dhvang))>0;        
                usign = sum(double([sign(dhbang)==1&sign(dhvang)==1,sign(dhbang)==-1&sign(dhvang)==-1])*[-1;1],2);
              case 2
                uquad = (sign(dhbang).*sign(dhvang))<0;            
                usign = sum(double([sign(dhbang)==1&sign(dhvang)==-1,sign(dhbang)==-1&sign(dhvang)==1])*[-1;1],2);            
            end            
            if e == 1,
                usign = ones(size(dhbang));                    
            end
            
            tferr = ferr{e}(ind & uquad).*usign(ind & uquad);            
            tdphz = dphz(ind & uquad);    
            thvangBinInd = hvangBinInd(ind & uquad);
            thbangBinInd = hbangBinInd(ind & uquad);
            tsmMask = smMask(ind & uquad);
            
            disp(['[Info] g: ',num2str(g),' e: ',num2str(e)]);
            for f = 1:edy,    
                for b = 1:edx
                    tic
                    for i = 1:numIter,
                        
                        indb = thvangBinInd==f & thbangBinInd==b & circshift(tsmMask, randsample(shifts,1));
                        indb(indb) = randn([sum(indb),1]) > -0.5;
                        pxavn(:,:,b,f,sts,e,g,i) =               ...
                            histcounts2(tferr(indb),          ...
                                        tdphz(indb),             ...
                                        ferrorBinEdges{e},       ...
                                        phzBins);
                    end
                    toc
                end
            end
        end
    end
end




npxav = bsxfun(@rdivide,pxavn,sum(pxavn));

qtls = [0.25,0.50,0.75];

qntlav = nan([3,size(npxav,2),size(npxav,3),size(npxav,4),size(npxav,5),size(npxav,6),size(npxav,7),size(npxav,8)]);
for p = 1:size(npxav,2),
    for i = 1:size(npxav,3),
        for s =1:size(npxav,4),
            disp(['[Info] sts: ',num2str(s)]);
            for t = 1:size(npxav,5),
                for f = 1:size(npxav,6),
                    for g= 1:size(npxav,7),
                        for e = 1:size(npxav,8),
                            try
                                [x, index] = unique(cumsum(npxav(:,p,i,s,t,f,g,e)));
                                qntlav(:,p,i,s,t,f,g,e) = interp1(x,ferrorBinCenters{f}(index),qtls); 
                            end
                        end
                    end
                end
            end
        end
    end
end




g = 1;
sts = 1;
e = 2
figure,
for f = 1:size(pxavn,4),    
    for b = 1:size(pxavn,3),    
        subplot2(size(pxavn,4),size(pxavn,3),f,b);
        imagesc(ferrorBinCenters{e},...
                [phzBinCenters;phzBinCenters+2*pi],...
                imgaussfilt(repmat(sq(mean(pxavn(:,:,b,f,sts,e,g,:),8))',[2,1]),[0.1,2]));
        axis('xy');
    end
end
ForAllSubplots('axis tight')
ForAllSubplots('xlim([-200,200]);');
ForAllSubplots('ylim([-pi,3*pi]);');
colormap jet




figure,
for f = 1:edy,    
    for b = 1:edx
        subplot2(5,5,f,b);
        plot(repmat(sq(qntlav(2,:,b,f,1,2,1,:)),[2,1]),[phzBinCenters;phzBinCenters+2*pi],'b')
        xlim([-100,100]);
        ylim([-pi,3*pi]);
    end
end

ForAllSubplots('Lines(0,[],''k'');');



figure,
for f = 1:edy,    
    for b = 1:edx
        subplot2(5,5,f,b);
        hold('on');
        plot(repmat(sq(mean(qntlav(3,:,b,f,:,2,1,:),8)),[2,1]),[phzBinCenters;phzBinCenters+2*pi])
        plot(repmat(sq(mean(qntlav(1,:,b,f,:,2,1,:),8)),[2,1]),[phzBinCenters;phzBinCenters+2*pi])        
        xlim([-100,200]);
        ylim([-pi,3*pi]);
    end
end
legend(stsLbls);

ForAllSubplots('Lines(0,[],''k'');');


quad = 2;
ledLoc = {'southwest','northeast'}; 
hba = 1:edx-1;

figure();
for p = [6,3],;
for hva = 1:4;
subplot(1,5,hva);
hold('on');
% qntlav( x, y, hba, hav, sts, fet, quad, itr)
plot( sq(mean(qntlav(2,p,hba,hva,:,1,quad,:),8,'omitnan')),sq(mean(qntlav(2,p,hba,hva,:,2,quad,:),8,'omitnan')))
xlim([-120,120]);
ylim([-120,120]);
legend(stsLbls,'Location',ledLoc{quad});
% $$$ plot(sq(mean(qntlav(2,6,hba,hva,:,1,1,:),8)),sq(mean(qntlav(2,6,hba,hva,:,2,1,:),8)))
% $$$ plot(sq(mean(qntlav(2,6,hba,hva,:,1,2,:),8)),sq(mean(qntlav(2,6,hba,hva,:,2,2,:),8)))
% $$$ plot(sq(mean(qntlav(2,3,hba,hva,:,1,2,:),8)),sq(mean(qntlav(2,3,hba,hva,:,2,2,:),8)))
end
end



hva = 1:4;
figure();
for hba = 1:8;
subplot(1,4,hba);
hold('on');
% pxav( x, y, hba, hav, sts, fet, quad, itr)
plot(sq(mean(qntlav(2,3,hba,hva,:,1,1,:),8)),sq(mean(qntlav(2,3,hba,hva,:,2,1,:),8)))
xlim([-60,120]);
ylim([-120,120]);
legend(stsLbls,'Location','southwest');
plot(sq(mean(qntlav(2,6,hba,hva,:,1,1,:),8)),sq(mean(qntlav(2,6,hba,hva,:,2,1,:),8)))
plot(sq(mean(qntlav(2,6,hba,hva,:,1,2,:),8)),sq(mean(qntlav(2,6,hba,hva,:,2,2,:),8)))
plot(sq(mean(qntlav(2,3,hba,hva,:,1,2,:),8)),sq(mean(qntlav(2,3,hba,hva,:,2,2,:),8)))
end


stt = 1:8:10000000;
ttt = [1:(10000000/8)]./250;
figure,
hold('on');
subplot(311);
hold('on');
plot(ttt,dhbang(stt));
plot(ttt,dhvang(stt).*10);
Lines([],0,'k');
subplot(312);
hold('on');
plot(ttt,dhrvf(stt));
plot(ttt,dhrvl(stt));
Lines([],0,'k');
subplot(313);
imagesc(ttt,1:8,dstcm(stt,:)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');




% pxav( x, y, hba, hav, sts, fet, quad, itr)
e = 2;
figure,
for p = 1:8,    
    for s = 1:5,
    subplot2(5,8,s,p);
    imagesc(sq(mean(qntlav(2,phzBinOrder(p),1:end-1,1:end,s,e,1,:),8))');
    if p==1,ylabel(stsLbls{s});end
    end
end
ForAllSubplots('caxis([-80,80]);');
ForAllSubplots('axis(''xy'')');
ForAllSubplots('colormap(''jet'')');
ForAllSubplots('Lines(0,[],''k'');');


indDD =   logical(dstcm(:,1))                             ...
            & any(dstcm(:,3),2)                  ...
            & dhdist<380; 



out = hist2([dhvang(ind),dhrvl(ind)],linspace(-0.2,0.2,50),linspace(-50,50,50));
figure();
imagescnan({linspace(-0.2,0.2,49),linspace(-50,50,49),abs(log10(out))'},[2,6]);

dhrvld = circshift(-dhrvl,30*8);
out = hist2([dhrvl(ind&dhbang<0.2&dhbang>-0.2),dhrvld(ind&dhbang<0.2&dhbang>-0.2)],linspace(-50,50,50),linspace(-50,50,50));
figure();
imagescnan({linspace(-50,50,49),linspace(-50,50,49),abs(log10(out))'},[2,6]);

figure,
for sts = 1:5,
ind =   logical(dstcm(:,1))                             ...
            & any(dstcm(:,stsGrps{sts}),2)                  ...
            & dpostI                                        ...
            & duincI ...
            & dhdist<380;

subplot(5,1,sts);
hist2([dhbang(ind),dhvang(ind)],linspace(-2,2,100),linspace(-0.1,0.1,100));
end

figure,
for sts = 1:5,
ind =   logical(dstcm(:,1))                             ...
            & any(dstcm(:,stsGrps{sts}),2)                  ...
            & dpostI                                        ...
            & duincI ...
            & dhdist<380;

subplot(5,1,sts);
hist2([log10(dxyvel(ind)),dhvang(ind)],linspace(-2,2,100),linspace(-0.1,0.1,100));
end

figure,
for sts = 1:5,
ind =   logical(dstcm(:,1))                             ...
            & any(dstcm(:,stsGrps{sts}),2)                  ...
            & dpostI                                        ...
            & duincI ...
            & dhdist<380;

subplot(5,1,sts);
hist2([log10(dxyvel(ind)),dhbang(ind)],linspace(-2,2,100),linspace(-pi/2,pi/2,100));
end





