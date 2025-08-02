% SECTIONS :
%     figure: HRVL -> Par-3 HBA
%     figure: HRVF -> Par-3 HBA
%
% TERMS :
%     
%    JPDF
global MTA_PROJECT_PATH
global MTA_PROJECT_REPORT_PATH

MTA_PROJECT_REPORT_PATH = '/storage/share/Projects/BehaviorPlaceCode';

create_directory(fullfile(MTA_PROJECT_REPORT_PATH,'egocentricPP'));

% Reverse ordered phase
phzOrder = [4,3,2,1,8,7,6,5];



%%%<<< testing effect of various behavioral variables on lateral phase precession

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
hbangBinInd = discretize(dhbang +0.4-0.8*double(~ismember(dtind,[3,4,5])),hbangBinEdges);

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


eds = 10;
hrvfBinEdges = linspace(-30,80,eds);
hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
hrvfBinInd = discretize(-dhrvf,hrvfBinEdges);

hrvlBinEdges = linspace(-60,60,eds);
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
hrvlBinInd = discretize(dhrvl,hrvlBinEdges);


% $$$ eds = 8;
% $$$ hrvfBinEdges = linspace(-30,80,eds);
% $$$ hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
% $$$ hrvfBinInd = discretize(-dhrvf,hrvfBinEdges);


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





eds = 10;
hrvfBinEdges = linspace(-30,80,eds);
hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
hrvfBinInd = discretize(-dfhrvf,hrvfBinEdges);

hrvlBinEdges = linspace(-60,60,eds);
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
hrvlBinInd = discretize(dfhrvl,hrvlBinEdges);

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
hvangBinEdges = [ -0.17,-0.018,0.018,0.17];
hvangBinEdges = [ -0.17,0.17];
%hvangBinEdges = linspace(-0.1,0.1,8);
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edy = numel(hvangBinCenters);
%hvangBinInd = discretize(circshift(dhvang,-25*8),hvangBinEdges);
hvangBinInd = discretize(dhvang,hvangBinEdges);
%hvangBinInd = discretize(circshift(dhvang,25*8),hvangBinEdges);

hbangBinEdges = linspace(-pi/3,pi/3,6);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.3-0.55*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);




ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<360;
%        & ~ismember(dtind,[3,4,5]);
%        & -dhrvf>=-10;




e = 2;
figure();
sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
tind = ind;% & ismember(hrollBinInd,5:7);
set(gcf,'Name','all');
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
ForAllSubplots('xlim([-250,250]);');
ForAllSubplots('colormap(''jet'');');
ForAllSubplots('Lines([],pi,''k'');');
ForAllSubplots('Lines(0,[],''k'');');
%ForAllSubplots('caxis(caxis()./1.5);');
%set(gcf,'Name','all');








hbangBinEdges = linspace(-1.2,1.2,10);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380; % ...


thrvfBinInd = hrvfBinInd(ind);
thrvlBinInd = hrvlBinInd(ind);
thvangBinInd = hvangBinInd(ind);

thbangBinInd = hbangBinInd(ind);
tdphz = dphz(ind);

%e = 1;
%figure();
%sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';

figure
sax = reshape(tight_subplot(8,6,0.001,0,0),[6,8])';
varlbls = {'hrvf','hrvl','hvang'};
for e = 1:2,
    tferr = ferr{e}(ind);
    for s = 1:numel(varlbls),
        switch varlbls{s},
          case 'hvang'
            var = 'hvang';
            hvangBinEdges = [-0.3,-0.06,-0.04,-0.024,-0.012,0.012,0.024,0.04,0.06,0.3];
            %hvangBinEdges = [-0.06,-0.04,-0.024,-0.018,-0.012,-0.006,0.006,0.012,0.018,0.024,0.04,0.06];
            hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
            edy = numel(hvangBinCenters);
            hvangBinInd = discretize(circshift(dhvang,0),hvangBinEdges);
            varBC = hvangBinCenters;
            thvangBinInd = hvangBinInd(ind);            
          case 'hrvf',
            var = 'hrvf';
            hrvfBinEdges = [-30,-15,-3,3,8,15,20,60,80];
            hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
            %edx = numel(hrvfBinCenters);
            edy = numel(hrvfBinCenters);
            hrvfBinInd = discretize(circshift(-dhrvf,0),hrvfBinEdges);
            varBC = hrvfBinCenters;
            thrvfBinInd = hrvfBinInd(ind);
          case 'hrvl'
            var = 'hrvl';
            hrvlBinEdges = [-60,-20,-15,-8,-3,3,8,15,20,60];
            hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
            %edx = numel(hrvlBinCennters);
            edy = numel(hrvlBinCenters);
            hrvlBinInd = discretize(circshift(-dhrvl,0),hrvlBinEdges);
            varBC = hrvlBinCenters;
            thrvlBinInd = hrvlBinInd(ind);
        end

        out = zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters),edx,edy]);
        for f = 1:edy,    
            for b = 1:edx
                %axes(sax(f,b));       
            switch var
              case 'hrvf'
                indb = thrvfBinInd==f &  thbangBinInd==b;
              case 'hrvl'
                indb = thrvlBinInd==f &  thbangBinInd==b;        
              case 'hvang'
                indb = thvangBinInd==f &  thbangBinInd==b;        
            end
            out(:,:,b,f) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
            end    
        end

        mout = zeros([1,numel(phzBinCenters),edx,edy]);
        for f = 1:edy,    
            for b = 1:edx
                for p = 1:numel(phzBinCenters),
                    try
                        [x, index] = unique(cumsum(out(:,p,b,f)./sum(out(:,p,b,f))),'first');
                        index(x==0|x==1) = [];
                        x(x==0|x==1) = [];
                        mout(:,p,b,f) = interp1(x,ferrorBinCenters{e}(index),0.5);
                    end
                end
            end;
        end

        for p = 1:8
            %subplot2(8,6,p,s+3*(e-1));
            axes(sax(p,s+3*(e-1)));
            imagesc(hbangBinCenters,varBC,sq(mout(:,phzOrder(p),:,:))');
            if e==1,
                caxis([-60,120]);        
            else,
                caxis([-60,60]);
            end
            colormap('jet');
        end
    end;
end;
figure();
imagesc(hbangBinCenters,varBC,sq(diff(mout(:,[3,7],:,:),1,2))');




%%%<<< HRVF and HRVL EPP over 3 partitions of HBA

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380; % ...

hvangBinEdges = [-0.3,-0.06,-0.04,-0.024,-0.012,0.012,0.024,0.04,0.06,0.3];
hvangBinEdges = [-0.06,-0.04,-0.024,-0.018,-0.012,-0.006,0.006,0.012,0.018,0.024,0.04,0.06];
hvangBinEdges = [-0.3,-0.06,-0.04,-0.024,-0.012,0.012,0.024,0.04,0.06,0.3];
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edy = numel(hvangBinCenters);
edx = numel(hvangBinCenters);
hvangBinInd = discretize(circshift(dhvang,0),hvangBinEdges);


hbangBinEdges = linspace(-1.2,1.2,10);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.45*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);


figure();
subplot(211);
hist(dhbang(~ismember(dtind,[3,4,5])&ind),linspace(-pi,pi,200));
circ_mean(dhbang(~ismember(dtind,[3,4,5])&ind))
subplot(212);
hist(dhbang(ismember(dtind,[3,4,5])&ind),linspace(-pi,pi,200));
circ_mean(dhbang(ismember(dtind,[3,4,5])&ind))

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


hbangBinEdges = linspace(-1.2,1.2,10);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.45*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);


var = 'hvang';
hvangBinEdges = [-0.3,-0.06,-0.04,-0.024,-0.012,0.012,0.024,0.04,0.06,0.3];
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edx = numel(hvangBinCenters);
hvangBinInd = discretize(circshift(dhvang,0),hvangBinEdges);
varBC = hvangBinCenters;
vx = hvangBinInd(ind);

var = 'hrvf';
hrvfBinEdges = [-30,-10,-4,-1,1,4,7,10,15,22,30,45,80];
hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
edx = numel(hrvfBinCenters);
edy = numel(hrvfBinCenters);
hrvfBinInd = discretize(circshift(-dfhrvf,0),hrvfBinEdges);
varBC = hrvfBinCenters;
thrvfBinInd = hrvfBinInd(ind);

var = 'hrvl';
hrvlBinEdges = [-60,-45,-30,-15,-10,-4,-1,1,4,10,15,30,45,60];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
%edx = numel(hrvlBinCennters);
edx = numel(hrvlBinCenters);
hrvlBinInd = discretize(circshift(dfhrvl,0),hrvlBinEdges);
varBC = hrvlBinCenters;
thrvlBinInd = hrvlBinInd(ind);


thbangBinInd = hbangBinInd(ind);
tdphz = dphz(ind);


figure,
subplot(211);
hist2([dhrvl(ind&~ismember(dtind,[3,4,5])&dhdist<300),-dhrvf(ind&~ismember(dtind,[3,4,5])&dhdist<300)],linspace(-100,100,100),linspace(-50,100,100))
caxis([0,1e3]);
subplot(212);
hist2([dhrvl(ind&ismember(dtind,[3,4,5])&dhdist<300),-dhrvf(ind&ismember(dtind,[3,4,5])&dhdist<300)],linspace(-100,100,100),linspace(-50,100,100))
caxis([0,1e3]);
shift = 0;
mout = zeros([1,numel(phzBinCenters),edx, edy, 3, 2]);
for e = 1:2;
    tferr = ferr{e}(ind);
    for g = 1:3,
        figure();
        sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
        for f = 1:edy,    for b = 1:edx
                axes(sax(f,b));               

                if g==1,
                    indb = thrvfBinInd==f & thrvlBinInd==b & ismember(thbangBinInd,1:4);
                    %indb = thrvfBinInd==f & thrvlBinInd==b & ismember(thbangBinInd,2:3);
                elseif g == 2,
                    indb = thrvfBinInd==f & thrvlBinInd==b & ismember(thbangBinInd,5);
                elseif g==3,
                    indb = thrvfBinInd==f & thrvlBinInd==b & ismember(thbangBinInd,6:9);
                end
                
                hold('on');
                out(:,:,b,f,g,e) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
                imagesc(ferrorBinEdges{e},    [phzBins,phzBins+2*pi],      imgaussfilt(repmat(out(:,:,b,f,g,e),[1,2]),[2,0.1])');
                axis('tight');        axis('xy');
% $$$         text(sax(f,b), -220, -pi,...
% $$$              ['hro:',num2str(round(hrollBinCenters(f),2)),...
% $$$               ' hva:',num2str(round(hvangBinCenters(b),2))],...
% $$$              'Color','w',...
% $$$              'Rotation',90);
% $$$         text(sax(f,b), -220, -pi,...
% $$$              ['hba:',num2str(round(hbangBinCenters(b),2)),...
% $$$               ' hro:',num2str(round(hrollBinCenters(f),2))],...
% $$$              'Color','w',...
% $$$              'Rotation',90);
                text(sax(f,b), -220, -pi,...
                     ['hvl:',num2str(round(hrvlBinCenters(b),2)),...
                      ' hvf:',num2str(round(hrvfBinCenters(f),2))],...
                     'Color','w',...
                     'Rotation',90);
% $$$                 text(sax(f,b), -220, -pi,...
% $$$                      ['hva:',num2str(round(hvangBinCenters(b),2)),...
% $$$                       ' hvf:',num2str(round(hrvfBinCenters(f),2))],...
% $$$                      'Color','w',...
% $$$                      'Rotation',90);
        end;end        

        ForAllSubplots('ylim([0,2*pi]);'); 
        ForAllSubplots('xlim([-250,250]);'); 
        ForAllSubplots('colormap(''jet'');');
        ForAllSubplots('Lines(0,[],''k'');');
        ForAllSubplots('Lines([],0,''k'');');  
        ForAllSubplots('Lines([],pi,''m'');');
        ForAllSubplots('Lines([],2*pi,''k'');');
    end
end

for f = 1:edy,    
    for b = 1:edx
        for p = 1:numel(phzBinCenters),
            for  g = 1:3,
                for e = 1:2,
                try
                    [x, index] = unique(cumsum(out(:,p,b,f,g,e)./sum(out(:,p,b,f,g,e))),'first');
                    index(x==0|x==1) = [];
                    x(x==0|x==1) = [];
                    mout(:,p,b,f,g,e) = interp1(x,ferrorBinCenters{e}(index),0.5);
                end
                end
            end
        end
    end;
end


figure();
for p = 1:8,
    subplot2(8,6,p,1);
    imagesc(hrvlBinCenters,hrvfBinCenters,sq(mout(:,phzOrder(p),:,:,1,1))');    colormap('jet');    caxis([-60,120]);
    subplot2(8,6,p,2);
    imagesc(hrvlBinCenters,hrvfBinCenters,sq(mout(:,phzOrder(p),:,:,2,1))');    colormap('jet');    caxis([-60,120]);    
    subplot2(8,6,p,3);
    imagesc(hrvlBinCenters,hrvfBinCenters,sq(mout(:,phzOrder(p),:,:,3,1))');    colormap('jet');    caxis([-60,120]);    
    subplot2(8,6,p,4);
    imagesc(hrvlBinCenters,hrvfBinCenters,sq(mout(:,phzOrder(p),:,:,1,2))');    colormap('jet');    caxis([-50,50]);
    subplot2(8,6,p,5);
    imagesc(hrvlBinCenters,hrvfBinCenters,sq(mout(:,phzOrder(p),:,:,2,2))');    colormap('jet');    caxis([-50,50]);    
    subplot2(8,6,p,6);
    imagesc(hrvlBinCenters,hrvfBinCenters,sq(mout(:,phzOrder(p),:,:,3,2))');    colormap('jet');    caxis([-50,50]);    
end

figure();
hold('on');
p = [':'];
subplot(311);
plot(reshape(mout(:,phzOrder(p),:,:,3,2),[],1),reshape(mout(:,phzOrder(p),:,:,3,1),[],1),'.c');
subplot(312);
plot(reshape(mout(:,phzOrder(p),:,:,2,2),[],1),reshape(mout(:,phzOrder(p),:,:,2,1),[],1),'.g');
subplot(313);
plot(reshape(mout(:,phzOrder(p),:,:,1,2),[],1),reshape(mout(:,phzOrder(p),:,:,1,1),[],1),'.m');
ForAllSubplots('xlim([-100,100])');
ForAllSubplots('ylim([-100,100])');

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


%%%>>>
%%%<<< HBA and HRVL EPP over 3 partitions of HRVF

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<320; % ...



hbangBinEdges = linspace(-1.2,1.2,10);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);
thbangBinInd = hbangBinInd(ind);

var = 'hrvf';
hrvfBinEdges = [-10,5,30,80];
hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
%edx = numel(hrvfBinCenters);
%edy = numel(hrvfBinCenters);
hrvfBinInd = discretize(circshift(-dfhrvf,0),hrvfBinEdges);
varBC = hrvfBinCenters;
thrvfBinInd = hrvfBinInd(ind);

var = 'hrvl';
hrvlBinEdges = [-60,-30,-15,-5,-1,1,5,15,30,60];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
%edx = numel(hrvlBinCennters);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(circshift(dfhrvl,0),hrvlBinEdges);
varBC = hrvlBinCenters;
thrvlBinInd = hrvlBinInd(ind);



tdphz = dphz(ind);
out = zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),edx, edy, 3, 2]);
for e = 1:2;
    tferr = ferr{e}(ind);
    for g = 1:3,
        figure();
        sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
        for f = 1:edy,    for b = 1:edx
                axes(sax(f,b));               

                if g==1,
                    indb = thrvlBinInd ==f & thbangBinInd==b & ismember(thrvfBinInd,1);
                elseif g == 2,
                    indb = thrvlBinInd ==f & thbangBinInd==b & ismember(thrvfBinInd,2);
                elseif g==3,
                    indb = thrvlBinInd ==f & thbangBinInd==b & ismember(thrvfBinInd,3);
                end
                
                hold('on');
                out(:,:,b,f,g,e) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
                imagesc(ferrorBinEdges{e}, [phzBins,phzBins+2*pi], imgaussfilt(repmat(out(:,:,b,f,g,e),[1,2]),[2,0.1])');
                axis('tight');        axis('xy');
                text(sax(f,b), -220, -pi,...
                     ['hba:',num2str(round(hbangBinCenters(b),2)),...
                      ' hvl:',num2str(round(hrvlBinCenters(f),2))],...
                     'Color','w',...
                     'Rotation',90);
        end;end        
        ForAllSubplots('ylim([0,2*pi]);'); 
        ForAllSubplots('xlim([-250,250]);'); 
        ForAllSubplots('colormap(''jet'');');
        ForAllSubplots('Lines(0,[],''k'');');
        ForAllSubplots('Lines([],0,''k'');');  
        ForAllSubplots('Lines([],pi,''m'');');
        ForAllSubplots('Lines([],2*pi,''k'');');
    end
end

mout = zeros([1,numel(phzBinCenters),edx, edy, 3, 2]);
for f = 1:edy,    
    for b = 1:edx
        for p = 1:numel(phzBinCenters),
            for  g = 1:3,
                for e = 1:2,
                try
                    [x, index] = unique(cumsum(out(:,p,b,f,g,e)./sum(out(:,p,b,f,g,e))),'first');
                    index(x==0|x==1) = [];
                    x(x==0|x==1) = [];
                    mout(:,p,b,f,g,e) = interp1(x,ferrorBinCenters{e}(index),0.5);
                end
                end
            end
        end
    end;
end


figure();
for p = 1:8,
    subplot2(8,6,p,1);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,1,1))');    colormap('jet');    caxis([-60,120]);
    subplot2(8,6,p,2);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,2,1))');    colormap('jet');    caxis([-60,120]);    
    subplot2(8,6,p,3);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,3,1))');    colormap('jet');    caxis([-60,120]);    
    subplot2(8,6,p,4);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,1,2))');    colormap('jet');    caxis([-50,50]);
    subplot2(8,6,p,5);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,2,2))');    colormap('jet');    caxis([-50,50]);    
    subplot2(8,6,p,6);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,3,2))');    colormap('jet');    caxis([-50,50]);    
end

figure();
hold('on');
p = [':'];
subplot(311);
plot(reshape(mout(:,phzOrder(p),:,:,3,2),[],1),reshape(mout(:,phzOrder(p),:,:,3,1),[],1),'.c');
subplot(312);
plot(reshape(mout(:,phzOrder(p),:,:,2,2),[],1),reshape(mout(:,phzOrder(p),:,:,2,1),[],1),'.g');
subplot(313);
plot(reshape(mout(:,phzOrder(p),:,:,1,2),[],1),reshape(mout(:,phzOrder(p),:,:,1,1),[],1),'.m');
ForAllSubplots('xlim([-100,100])');
ForAllSubplots('ylim([-100,100])');

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


%%%>>>

%%%<<< HBA and HRVL EPP

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<320; % ...

var = 'hbang';
hbangBinEdges = linspace(-1.2,1.2,8);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);
thbangBinInd = hbangBinInd(ind);

var = 'hrvl';
hrvlBinEdges = [-60,-45,-30,-20,-10,-5,-1,1,5,10,20,30,45,60];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
%edx = numel(hrvlBinCennters);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(circshift(dfhrvl,0),hrvlBinEdges);
varBC = hrvlBinCenters;
thrvlBinInd = hrvlBinInd(ind);

var = 'hrvf';
hrvlBinEdges = [-5,5,10,15,20,30,45,60,80];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(circshift(dfhrvf,0),hrvlBinEdges);
varBC = hrvlBinCenters;
thrvlBinInd = hrvlBinInd(ind);



tdphz = dphz(ind);


out =  zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters),edx, edy, 2]);
for e = 1:4;
    tferr = ferr{e}(ind);
    figure();
    sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
    for f = 1:edy,    for b = 1:edx
            axes(sax(f,b));               
            hold('on');
            indb = thrvlBinInd==f & thbangBinInd==b;
            out(:,:,b,f,e) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
            imagesc(ferrorBinEdges{e},    [phzBins,phzBins+2*pi],      imgaussfilt(repmat(out(:,:,b,f,e),[1,2]),[5,0.1])');
            axis('tight');        axis('xy');
            text(sax(f,b), -220, -pi,...
                 ['hba:',num2str(round(hbangBinCenters(b),2)),...
                  ' hvl:',num2str(round(hrvlBinCenters(f),2))],...
                 'Color','w',...
                 'Rotation',90);
    end;end        
    ForAllSubplots('ylim([0,2*pi]);'); 
    ForAllSubplots('xlim([-250,250]);'); 
    ForAllSubplots('colormap(''jet'');');
    ForAllSubplots('Lines(0,[],''k'');');
    ForAllSubplots('Lines([],0,''k'');');  
    ForAllSubplots('Lines([],pi,''m'');');
    ForAllSubplots('Lines([],2*pi,''k'');');
end

clear('mout');
mout = [];
for f = 1:edy,    
    for b = 1:edx
        for e = 1:2,
            mout(:,:,b,f,e) = imgaussfilt(repmat(out(:,:,b,f,e),[1,2]),[20,0.1]);
        end
    end
end
mout = mout(:,[9:12,5:8],:,:,:);
[~,mout] = max(mout);
mout = ferrorBinCenters{1}(mout);

mout = sum(bsxfun(@rdivide,out,sum(out)).*repmat(ferrorBinCenters{1}',[1,numel(phzBinCenters),edx,edy,2]));

mout = zeros([1,numel(phzBinCenters),edx, edy, 2]);
for f = 1:edy,    
    for b = 1:edx
        for p = 1:numel(phzBinCenters),
            for e = 1:2,
                try
                    [x, index] = unique(cumsum(out(:,p,b,f,e)./sum(out(:,p,b,f,e))),'first');
                    index(x==0|x==1) = [];
                    x(x==0|x==1) = [];
                    mout(:,p,b,f,e) = interp1(x,ferrorBinCenters{e}(index),0.5);
                end
            end
        end
    end
end


figure();
for p = 1:8,
    subplot2(8,2,p,1);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,1))');    colormap('jet');    caxis([-60,120]);
    subplot2(8,2,p,2);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,2))');    colormap('jet');    caxis([-60,60]);
end

figure();
hold('on');
p = [':'];
subplot(311);
plot(reshape(mout(:,phzOrder(p),:,:,3,2),[],1),reshape(mout(:,phzOrder(p),:,:,3,1),[],1),'.c');
subplot(312);
plot(reshape(mout(:,phzOrder(p),:,:,2,2),[],1),reshape(mout(:,phzOrder(p),:,:,2,1),[],1),'.g');
subplot(313);
plot(reshape(mout(:,phzOrder(p),:,:,1,2),[],1),reshape(mout(:,phzOrder(p),:,:,1,1),[],1),'.m');
ForAllSubplots('xlim([-100,100])');
ForAllSubplots('ylim([-100,100])');

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


%%%>>>




%%%<<< HBA and HRVL EPP

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380; % ...

var = 'dfet';
hbangBinEdges = linspace(-1.2,1.2,10);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);
thbangBinInd = hbangBinInd(ind);

var = 'fhrvl';
hrvlBinEdges = [-60,-45,-30,-20,-10,-5,-1,1,5,10,20,30,45,60];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
%edx = numel(hrvlBinCennters);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(circshift(dfhrvl,0),hrvlBinEdges);
%hrvlBinInd = discretize(circshift(dhrvl,0),hrvlBinEdges);
varBC = hrvlBinCenters;
thrvlBinInd = hrvlBinInd(ind);



hbangBinEdges = linspace(-1.5,0.2,10);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(dfet,hbangBinEdges);
thbangBinInd = hbangBinInd(ind);

hrvlBinEdges = [-30,-10,-4,-1,1,4,7,10,15,22,30,45,80];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(-dfhrvf,hrvlBinEdges);
thrvlBinInd = hrvlBinInd(ind);



tdphz = dphz(ind);


out =  zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters),edx, edy, 2]);
for e = 1:2;
    tferr = ferr{e}(ind);
    figure();
    sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
    for f = 1:edy,    for b = 1:edx
            axes(sax(f,b));               
            hold('on');
            indb = thrvlBinInd==f & thbangBinInd==b;
            out(:,:,b,f,e) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
            imagesc(ferrorBinEdges{e},    [phzBins,phzBins+2*pi],      imgaussfilt(repmat(out(:,:,b,f,e),[1,2]),[5,0.1])');
            axis('tight');        axis('xy');
            text(sax(f,b), -220, -pi,...
                 ['hba:',num2str(round(hbangBinCenters(b),2)),...
                  ' hvl:',num2str(round(hrvlBinCenters(f),2))],...
                 'Color','w',...
                 'Rotation',90);
    end;end        
    ForAllSubplots('ylim([0,2*pi]);'); 
    ForAllSubplots('xlim([-250,250]);'); 
    ForAllSubplots('colormap(''jet'');');
    ForAllSubplots('Lines(0,[],''k'');');
    ForAllSubplots('Lines([],0,''k'');');  
    ForAllSubplots('Lines([],pi,''m'');');
    ForAllSubplots('Lines([],2*pi,''k'');');
end

% $$$ clear('mout');
% $$$ mout = [];
% $$$ for f = 1:edy,    
% $$$     for b = 1:edx
% $$$         for e = 1:2,
% $$$             mout(:,:,b,f,e) = imgaussfilt(repmat(out(:,:,b,f,e),[1,2]),[20,0.1]);
% $$$         end
% $$$     end
% $$$ end
% $$$ mout = mout(:,[9:12,5:8],:,:,:);
% $$$ [~,mout] = max(mout);
% $$$ mout = ferrorBinCenters{1}(mout);
% $$$ mout = sum(bsxfun(@rdivide,out,sum(out)).*repmat(ferrorBinCenters{1}',[1,numel(phzBinCenters),edx,edy,2]));

mout = zeros([1,numel(phzBinCenters),edx, edy, 2]);
for f = 1:edy,    
    for b = 1:edx
        for p = 1:numel(phzBinCenters),
            for e = 1:2,
                try
                    [x, index] = unique(cumsum(out(:,p,b,f,e)./sum(out(:,p,b,f,e))),'first');
                    index(x==0|x==1) = [];
                    x(x==0|x==1) = [];
                    mout(:,p,b,f,e) = interp1(x,ferrorBinCenters{e}(index),0.5);
                end
            end
        end
    end
end


figure();
for p = 1:8,
    subplot2(8,2,p,1);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,1))');    colormap('jet');    caxis([-60,120]);
    subplot2(8,2,p,2);
    imagesc(hbangBinCenters,hrvlBinCenters,sq(mout(:,phzOrder(p),:,:,2))');    colormap('jet');    caxis([-60,60]);
end

figure();
hold('on');
p = [':'];
subplot(311);
plot(reshape(mout(:,phzOrder(p),:,:,3,2),[],1),reshape(mout(:,phzOrder(p),:,:,3,1),[],1),'.c');
subplot(312);
plot(reshape(mout(:,phzOrder(p),:,:,2,2),[],1),reshape(mout(:,phzOrder(p),:,:,2,1),[],1),'.g');
subplot(313);
plot(reshape(mout(:,phzOrder(p),:,:,1,2),[],1),reshape(mout(:,phzOrder(p),:,:,1,1),[],1),'.m');
ForAllSubplots('xlim([-100,100])');
ForAllSubplots('ylim([-100,100])');

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


%%%>>>


%%%<<< HRVL -> Par-3 HBA
%hbangBinEdges = linspace(-1.2,1.2,10);
hbangBinEdges = [-1.2,-0.15,0.15,1.2];
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
hbangBinInd = discretize(-(dhbang+0.2-0.45*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);
thbangBinInd = hbangBinInd(ind);

var = 'hrvl';
hrvlBinEdges = [-60,-45,-30,-20,-15,-10,-5,-1,1,5,10,15,20,30,45,60];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edx = numel(hrvlBinCenters);
hrvlBinInd = discretize(dfhrvl,hrvlBinEdges);
varBC = hrvlBinCenters;
thrvlBinInd = hrvlBinInd(ind);

tphz = dphz(tind);

out = zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters),edx,  3, 2]);
hbangPar = {1:4,5,6:9};
hbangPar = {1,2,3};
tic
for e = 1:2;
    tferr = ferr{e}(ind);
    for g = 1:3,
        for b = 1:edx
            indb = thrvlBinInd==b & ismember(thbangBinInd,hbangPar{g});
            out(:,:,b,g,e) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
        end
    end
end
toc

mout = zeros([numel(phzBinCenters),edx,  3, 2]);
for b = 1:edx
    for p = 1:numel(phzBinCenters),
        for g = 1:3,
            for e = 1:2,
                try
                    [x, index] = unique(cumsum(out(:,p,b,g,e)./sum(out(:,p,b,g,e))),'first');
                    index(x==0|x==1) = [];
                    x(x==0|x==1) = [];
                    mout(p,b,g,e) = interp1(x,ferrorBinCenters{e}(index),0.5);
                end
            end
        end
    end
end;

figure();
for g = 1:3,
    for e = 1:2,
        subplot2(2,3,e,g);
        imagesc(hrvlBinCenters,phzBinCenters(phzOrder)+2*pi*double(phzBinCenters(phzOrder)<0),mout(phzOrder,:,g,e));
        colormap('jet');
        axis('xy');
        if e==1,
            caxis([-60,120]);
        else
            caxis([-50,50]);
        end
    end
end


%%%>>>

%%%<<< HRVF -> Par-3 HBA
hbangBinEdges = linspace(-1.2,1.2,10);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.45*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);

var = 'hrvf';
hrvfBinEdges = [-30,-10,-4,-1,1,4,7,10,15,22,30,45,80];
hrvfBinCenters = mean([hrvfBinEdges(2:end); hrvfBinEdges(1:end-1)]);
edy = numel(hrvfBinCenters);
hrvfBinInd = discretize(circshift(-dhrvf,0),hrvfBinEdges);
varBC = hrvfBinCenters;
thrvfBinInd = hrvfBinInd(ind);

out = zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters),edx,  3, 2]);
hbangPar = {1:4,5,6:9};
tic
for e = 1:2;
    tferr = ferr{e}(ind);
    for g = 1:3,
        for b = 1:edx
            indb = thrvfBinInd==b & ismember(thbangBinInd,hbangPar{g});
            out(:,:,b,g,e) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
        end
    end
end
toc

mout = zeros([numel(phzBinCenters),edx,  3, 2]);
for b = 1:edx
    for p = 1:numel(phzBinCenters),
        for  g = 1:3,
            for e = 1:2,
                try
                    [x, index] = unique(cumsum(out(:,p,b,g,e)./sum(out(:,p,b,g,e))),'first');
                    index(x==0|x==1) = [];
                    x(x==0|x==1) = [];
                    mout(p,b,g,e) = interp1(x,ferrorBinCenters{e}(index),0.5);
                end
            end
        end
    end
end;


figure();
for g = 1:3,
    for e = 1:2,
        subplot2(2,3,e,g);
        imagesc(hrvfBinCenters,phzBinCenters(phzOrder)+2*pi*double(phzBinCenters(phzOrder)<0),mout(phzOrder,:,g,e));
        colormap('jet');
        axis('xy');
        if e==1,
            caxis([-60,120]);
        else
            caxis([-50,50]);
        end
    end
end


%%%>>>


%%%<<< something which needs a name later

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)              ...
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

%%%>>>

%%%<<< figure part hrvl and hba
ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[10])),2)            ...
        & duincI                                        ...
        & dpostI                                        ...        
        & dhdist<380;


var = 'hbang';
%hbangBinEdges = [-1.2,-0.6,-0.2,0.2,0.6,1.2];
%hbangBinEdges = linspace(-1.2,1.2,12);
hbangBinEdges = linspace(-0.5,0.5,4);
%hbangBinEdges = [-1.2,1.2];
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);
thbangBinInd = hbangBinInd(ind);
%tdhbang = -(dhbang+0.2-0.45*double(~ismember(dtind,[3,4,5])));
%tdhbang = tdhbang(ind);

var = 'hrvl';
hrvlBinEdges = [-60,-15,-5,5,15,60];    
%hrvlBinEdges = [-60,-5,5,60];    
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
%edx = numel(hrvlBinCennters);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(dfhrvl,hrvlBinEdges);
thrvlBinInd = hrvlBinInd(ind);
%tdfhrvl = dfhrvl(ind);

var = 'fhrvf';
hrvlBinEdges = [-10,0,5,15,30,100];
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(-dfhrvf,hrvlBinEdges);
thrvlBinInd = hrvlBinInd(ind);


var = 'hbpch';
hrvlBinEdges = linspace(-1.5,0.5,26);
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(dfet,hrvlBinEdges);
thrvlBinInd = hrvlBinInd(ind);


var = 'hpch';
hrvlBinEdges = linspace(-1.5,0.5,16);
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(dpch,hrvlBinEdges);
thrvlBinInd = hrvlBinInd(ind);

var = 'bpch';
hrvlBinEdges = linspace(-.1,0.4,16);
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(dbch,hrvlBinEdges);
thrvlBinInd = hrvlBinInd(ind);

var = 'rhm';
rhmBinEdges = phzBins;
rhmBinCenters = mean([rhmBinEdges(2:end); rhmBinEdges(1:end-1)]);
edy = numel(rhmBinCenters);
rhmBinInd = discretize(circshift(drhm,0),rhmBinEdges);
trhmBinInd = rhmBinInd(ind);


figure,hist2([-(dhbang(ind)+0.2-0.4*double(~ismember(dtind(ind),[3,4,5]))),drhm(ind)],linspace(-1.2,1.2,100),linspace(-pi,pi,16));

hrvlBinEdges = linspace(30,110,26);
hrvlBinCenters = mean([hrvlBinEdges(2:end); hrvlBinEdges(1:end-1)]);
edy = numel(hrvlBinCenters);
hrvlBinInd = discretize(dhz,hrvlBinEdges);
thrvlBinInd = hrvlBinInd(ind);


tdphz = dphz(ind);
tdfhrvf= -dfhrvf(ind);
tsmMask = smMask(ind);
% compute JPDF(FERROR, PHZ | HBA, HRVL)
out =  zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters),edx, edy, 2,10]);
for e = 1:2;
    tferr = ferr{e}(ind);
    for f = 1:edy,
        for b = 1:edx
            for i = 1:10
            indb = trhmBinInd==f & thbangBinInd==b & tsmMask;
            indb(randn([sum(indb),1]) > -0.7) = false;
            out(:,:,b,f,e,i) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
            end
        end        
    end
end

% plot JPDF(FERROR, PHZ | HBA, HRVL)
for e = 1:2;
    figure();
    sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
    for f = 1:edy,    for b = 1:edx
            axes(sax(f,b));               
            hold('on');
            imagesc(ferrorBinEdges{e},                                           ...
                    [phzBins,phzBins+2*pi],                                      ...
                    imgaussfilt(repmat(mean(out(:,:,b,f,e,:),ndims(out)),[1,2]),[2,0.1])');
            axis('tight');        axis('xy');
            text(sax(f,b), -220, -pi,...
                 ['hba:',num2str(round(hbangBinCenters(b),2)),...
                  ' hvl:',num2str(round(rhmBinCenters(f),2))],...
                 'Color','w',...
                 'Rotation',90);
    end;end        
    ForAllSubplots('ylim([0,2*pi]);'); 
    if e==1,
        ForAllSubplots('xlim([-200,350]);'); 
    else,
        ForAllSubplots('xlim([-300,300]);');         
    end
    ForAllSubplots('colormap(''jet'');');
    ForAllSubplots('Lines(0,[],''k'');');
    ForAllSubplots('Lines([],0,''k'');');  
    ForAllSubplots('Lines([],pi,''m'');');
    ForAllSubplots('Lines([],2*pi,''k'');');
end

mout = zeros([1,numel(phzBinCenters),edx, edy, 2]);
for f = 1:edy,    
    for b = 1:edx
        for p = 1:numel(phzBinCenters),
            for e = 1:2,
                for i = 1:10,
                try
                    [x, index] = unique(cumsum(out(:,p,b,f,e,i)./sum(out(:,p,b,f,e,i))),'last');
                    index(x==0|x==1) = [];
                    x(x==0|x==1) = [];
                    mout(:,p,b,f,e,i) = interp1(x,ferrorBinCenters{e}(index),0.5);
                end
                end
            end
        end
    end
end
 

figure();
for p = 1:8,
    subplot2(8,2,p,1);
    imagesc(hbangBinCenters,rhmBinCenters,sq(mean(mout(:,phzOrder(p),:,:,1,:),ndims(mout)))');  colormap('jet');  caxis([-40,100]);  axis('xy');
    subplot2(8,2,p,2);
    imagesc(hbangBinCenters,rhmBinCenters,sq(mean(mout(:,phzOrder(p),:,:,2,:),ndims(mout)))');  colormap('jet');  caxis([-50,50]);  axis('xy');
end 

figure,plot(sq(mout(1,phzOrder(:),1,:,1))')

cmap = hsv(15);
figure();
hold('on');
for p = 2:14,
plot(sq(mout(1,:,1,p,1)),'Color',cmap(p,:));
end


figure();imagesc(phzBinCenters(phzOrder),hrvlBinCenters,sq(mout(1,phzOrder,1,:,1))')
axis('xy');


sla = cell([1,2]);
%[sla{:}] = meshgrid([-100,0,100],[-100,0,100]);
[sla{:}] = meshgrid([-150,-75,0,75,150],[-150,-75,0,75,150]);
sla = cf(@(s) reshape(s,[],1), sla);


cmap = cool(7);
figure();
hold('on');
for p = 6,
    quiver(sla{[2,1]},...
           reshape(mout(:,phzOrder(p),:,:,2),[],1),...
           reshape(mout(:,phzOrder(p),:,:,1),[],1),...
           0,...
           'Color',cmap(p,:));
end
af(@(a) set(a,'ShowArrowHead','off'), findobj(gcf,'Type','quiver'));
set(gca,'YTick',[-100,0,100]);
set(gca,'XTick',[-100,0,100]);
grid(gca,'on');
daspect([1,1,1])


cmap = cool(7);
figure();
hold('on');
for p = 2,
    quiver(sla{[2,1]},                                       ...
           reshape(mout(:,phzOrder(p),:,:,2),[],1)./1.7,     ...
           reshape(mout(:,phzOrder(p),:,:,1),[],1)./1.7,     ...
           0,                                                ...
           'Color',cmap(p,:));
end
af(@(a) set(a,'ShowArrowHead','off'), findobj(gcf,'Type','quiver'));
set(gca,'YTick',[-75,0,75]);
set(gca,'XTick',[-75,0,75]);


cmap = cool(5);
figure,
hold('on');
quiver(sla{[2,1]},reshape(diff(mout(:,[6,8],:,:,2),1,2),[],1),reshape(diff(mout(:,[6,8],:,:,1),1,2),[],1),'Color',cmap(1,:));
quiver(sla{[2,1]},reshape(diff(mout(:,[6,1],:,:,2),1,2),[],1),reshape(diff(mout(:,[6,1],:,:,1),1,2),[],1),'Color',cmap(2,:));
quiver(sla{[2,1]},reshape(diff(mout(:,[6,2],:,:,2),1,2),[],1),reshape(diff(mout(:,[6,2],:,:,1),1,2),[],1),'Color',cmap(3,:));
quiver(sla{[2,1]},reshape(diff(mout(:,[6,3],:,:,2),1,2),[],1),reshape(diff(mout(:,[6,3],:,:,1),1,2),[],1),'Color',cmap(4,:));
quiver(sla{[2,1]},reshape(diff(mout(:,[6,4],:,:,2),1,2),[],1),reshape(diff(mout(:,[6,4],:,:,1),1,2),[],1),'Color',cmap(5,:));
af(@(a) set(a,'ShowArrowHead','off'), findobj(gcf,'Type','quiver'));
%%%>>>



%%%<<< figure part hrvl and hba mirrored

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)            ...
        & duincI                                        ...
        & dpostI                                        ...        
        & dhdist<380;

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[5])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[9])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;


ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[10])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;


ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[4,6])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;
ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[6])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;


ind =   logical(dstcm(:,1))                             ...
        & logical(dstcm(:,[5])) & logical(dstcm(:,[10])) ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,5])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;


ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[4,6])),2)              ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;


cdhroll = (dhroll-0.26*double(~ismember(dtind,[3,4,5])));

figure,
subplot(311);
hist(cdhroll(dhrvf>5&logical(dstcm(:,5))&ismember(dtind,[3,4,5])),1000)
subplot(312);
hist(cdhroll(dhrvf>5&logical(dstcm(:,5))&~ismember(dtind,[3,4,5])),1000)
subplot(313);
hist(cdhroll(dhrvf>5&logical(dstcm(:,5))),1000)
Lines(mean(cdhroll(dhrvf>5&logical(dstcm(:,5)))),[],'r');

chbang = -(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5])));


figure();
hold('on');
plot(dct.roll{1});
plot(dct.hbang{1}.data)
plot(dct.fhrvfl{1}(:,2)/50)


t = 14950;
t = 62375;
t = 4855;
t = 2260;
t = 2402;
mars = {'spine_lower','pelvis_root','spine_middle','head_back','head_left','hcom','head_right','nose'};
figure();
hold('on');
plot3(dct.xyz{1}(t,mars,1),dct.xyz{1}(t,mars,2),dct.xyz{1}(t,mars,3),'.-');
%plot3(dct.xyz{1}(t,'head_front',1),dct.xyz{1}(t,'head_front',2),dct.xyz{1}(t,'head_front',3),'.r','MarkerSize',20);
plot3(dct.xyz{1}(t,'nose',1),dct.xyz{1}(t,'nose',2),dct.xyz{1}(t,'nose',3),'.m','MarkerSize',20);
plot3(dct.xyz{1}(t,'head_back',1),dct.xyz{1}(t,'head_back',2),dct.xyz{1}(t,'head_back',3),'.c','MarkerSize',20);
daspect([1,1,1])


% cdhroll: neg -> left; pos -> right
% chbang:  neg -> left; pos -> right
% hvl:     neg -> left; pos -> right

%figure();hist2([chbang(ind),cdhroll(ind)],linspace(-1.2,1.2,50),linspace(-0.6,0.6,50));

%%%<<< COMPUTE EPP for up to 3 behavioral variables







ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)            ...
        & duincI                                        ...
        & dpostI                                        ...        
        & dhdist<380;


chroll = (dhroll-0.26*double(~ismember(dtind,[3,4,5])));
chbang = -(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5])));
chvang = dhvang;
cfhrvl = dfhrvl;
lfErr = ferr([1:2]);
v = cfhrvl;
a = chbang;
chbang  ( a<0 ) = -chbang   ( a<0 );
lfErr{2}( a<0 ) = -lfErr{2} ( a<0 );
cfhrvl  ( a<0 ) = -cfhrvl   ( a<0 );
chroll  ( a<0 ) = -chroll  ( a<0 );
chvang  ( a<0 ) = -chvang   ( a<0 );
% $$$ cfhrvl( a<0 & v>0) = -cfhrvl( a<0 & v>0);
% $$$ cfhrvl( a<0 & v<0) = -cfhrvl( a<0 & v<0);



vlb = {};
vbe = {};
vbc = {};
vbi = {};

vlb{end+1} = 'hbang';
vbe{end+1} = [0,0.12,0.26,0.44,0.7,1.2]; %v1
%vbe{end+1} = [0,0.055,0.123,0.18,0.25,0.325,0.4,0.49,0.613,0.76,0.97,1.2]; %v2
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);

% $$$ vlb{end+1} = 'hvang';
% $$$ %vbe{end+1} = [-0.2,-0.06,-0.012,0.012,0.06,0.2]; %v1
% $$$ vbe{end+1} = [-0.25,-0.15,-0.1,-0.07,-0.04,-0.015,0.015,0.04,0.07,0.1,0.15,0.25]; % v2
% $$$ vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
% $$$ vbi{end+1} = discretize(chvang(ind),vbe{end});

vlb{end+1} = 'hroll';
vbe{end+1} = linspace(-0.5,0.5,6); % v1
%vbe{end+1} = [-0.6,-0.3,-0.2,-0.15,-0.075,-0.02,0.02,0.075,0.15,0.2,0.3,0.6];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);

vlb{end+1} = 'hrvl';
vbe{end+1} = [-80,-10,-3,3,10,80]; %v1
%vbe{end+1} = [-60,-20,-10,-6,-3,-1,1,3,6,10,20,60]; % v2
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);

vlb{end+1} = 'fhrvf';
vbe{end+1} = [-3,3,10,20,40,80]; % v1
%vbe{end+1} = [-20,-6,-3,-1,1,4,8,12,18,28,45,80]; % v2
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);

%vrs = 'v1';
%vrs = 'v2';
vrs = 'v3';

%grps = {[1,2,3],[1,2,4],[1,2,5],[1,3,4],[1,3,5],[1,4,5],[2,3,4],[2,3,5],[2,4,5],[3,4,5]}; % v1
%grps = {[1,3,4],[1,3,5],[1,4,5]}; % v2
grps = {[1,2,3],[1,2,4],[1,3,4]};
stgrps = {[3:6],4,6,3,5,9,10};

nbins = 5;
mJpdf =  zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),nbins, nbins, nbins, 2,numel(grps),numel(stgrps),100]);
smplCnt =  zeros([nbins, nbins, nbins, 2,numel(grps),numel(stgrps),100]);
mmJpdf = zeros([3,numel(phzBinCenters),nbins, nbins, nbins, 2,numel(grps),numel(stgrps),100]);



%tdphz = dphz(ind);
  


% $$$ figure();
% $$$ stgrps = [4,6,3,5,9,10];
% $$$ nbins = 11;
% $$$ ny = 4;
% $$$ for s = 1:numel(stgrps)
% $$$ j = 1;
% $$$ ind =   logical(dstcm(:,1))                             ...
% $$$         & dstcm(:,stgrps(s))==stgrps(s)                 ...
% $$$         & duincI                                        ...
% $$$         & dpostI                                        ...        
% $$$         & dhdist<380;
% $$$ subplot2(numel(stgrps),ny,s,j);j = j+1;
% $$$ hist2([chbang(ind),chroll(ind)],linspace(0,1.2,nbins),linspace(-0.6,0.6,nbins));
% $$$ title('hba x hrl');
% $$$ ylabel(states{stgrps(s)})
% $$$ subplot2(numel(stgrps),ny,s,j);j = j+1;
% $$$ hist2([chbang(ind),cfhrvl(ind)],linspace(0,1.2,nbins),linspace(-60,60,nbins));
% $$$ title('hba x hvl');
% $$$ subplot2(numel(stgrps),ny,s,j);j = j+1;
% $$$ hist2([chbang(ind),dfhrvf(ind)],linspace(0,1.2,nbins),linspace(-20,80,nbins));
% $$$ title('hba x hvf');
% $$$ subplot2(numel(stgrps),ny,s,j);j = j+1;
% $$$ hist2([dfhrvf(ind),cfhrvl(ind)],linspace(-20,80,nbins),linspace(-60,60,nbins));
% $$$ title('hvf x hvl');
% $$$ end




for s = 1:numel(stgrps),

    
    ind =   logical(dstcm(:,1))                             ...
            & any(logical(dstcm(:,stgrps{s})),2)            ...
            & duincI                                        ...
            & dpostI                                        ...        
            & dhdist<380;

    vbi{1} = discretize(chbang(ind),vbe{1});
    vbi{2} = discretize(chroll(ind),vbe{2});
    vbi{3} = discretize(cfhrvl(ind),vbe{3});
    vbi{4} = discretize(dfhrvf(ind),vbe{4});



    lfErrInd = cf(@(e) discretize(e,ferrorBinEdges{1}), lfErr);
    phzBinInds = discretize(dphz,phzBins);
    tdphz = phzBinInds(ind);
    numPhzBins = numel(phzBinCenters);
    numFErBins = numel(ferrorBinCenters{1});

    for g = 1:numel(grps),

        xi = vbi{grps{g}(1)};
        xb = vbc{grps{g}(1)};
        xc = numel(xb);
        yi = vbi{grps{g}(2)};
        yb = vbc{grps{g}(2)};
        yc = numel(yb);
        zi = vbi{grps{g}(3)};
        zb = vbc{grps{g}(3)};
        zc = numel(zb);

        shifts = 0:8:2^8;
        tsmMask = smMask(ind);
        for e = 1:2;
            tferr = lfErrInd{e}(ind);
            tind = ~isnan(tferr+tdphz);
            tind(tind) = logical(tferr(tind)+tdphz(tind));

            % compute JPDF(FERROR, PHZ | HBA, HRVL)
            for x = 1:xc,
                for y = 1:yc,
                    for z = 1:zc,
                        disp(['s: ', num2str(s),'  g: ',num2str(g),'   x: ',num2str(x),'  y: ',num2str(y),'  z: ',num2str(z),'  e: ',num2str(e)]);
                        tic
                        for i = 1:100,

                            indb =     (xi==x) ...
                                     & (yi==y) ...
                                     & (zi==z) ...                                  
                                     & circshift(tsmMask, randsample(shifts,1)) ...
                                     & tind;
% $$$                             indb =   (xi==x-1|xi==x|xi==x+1) ...
% $$$                                      & (yi==y-1|yi==y|yi==y+1) ...
% $$$                                      & (zi==z-1|zi==z|zi==z+1) ...                                  
% $$$                                      & circshift(tsmMask, randsample(shifts,1)) ...
% $$$                                      & tind;
                            indb(randn([sum(indb),1]) > 0) = false;
                            smplCnt(x,y,z,e,g,s,i) = sum(indb);                        
                            mJpdf(:,:,x,y,z,e,g,s,i) = full(sparse(tferr(indb), tdphz(indb), 1, numFErBins, numPhzBins));
% $$$                         for p = 1:8,
% $$$                             for f = 1:99,
% $$$                                 mJpdf(:,:,x,y,z,e,g,i) = hist2([tferr(indb),tdphz(indb)],...
% $$$                                                      ferrorBinEdges{e},...
% $$$                                                      phzBins);
                        end
                        toc
                    end
                end
            end
        end
    end
end

%%%<<< Compute mmJpdf
for s = 1:size(mJpdf,8)
    for g = 1:size(mJpdf,7)
        for x = 1:size(mJpdf,3)
            tic
            for y = 1:size(mJpdf,4),
                for z = 1:size(mJpdf,5),
                    for p = 1:size(mJpdf,2),
                        for e = 1:size(mJpdf,6),
                            for i = 1:size(mJpdf,9)
                                try
                                    [c, index] = unique(cumsum(mJpdf(:,p,x,y,z,e,g,s,i)./sum(mJpdf(:,p,x,y,z,e,g,s,i))),'last');
                                    index(c==0|c==1) = [];
                                    c(c==0|c==1) = [];
                                    mmJpdf(:,p,x,y,z,e,g,s,i) = interp1(c,ferrorBinCenters{e}(index),[0.25,0.5,0.75]);
                                end
                            end
                        end
                    end
                end
            end
            toc
        end
    end
end
%%%>>>


save(fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_',vrs,'.mat']),...
     'mmJpdf',                                                                           ...
     'mJpdf',                                                                            ...
     'vlb',                                                                              ...
     'vbe',                                                                              ...
     'vbc',                                                                              ...
     'phzBins',                                                                          ...
     'grps',                                                                             ...
     '-v7.3'                                                                             ...
);
load(fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_',vrs,'.mat']));


%%%>>>

%labels    'hbang'    'hvang'    'hroll'    'hrvl'    'fhrvf'

grps = {[1,2,3], ... hba, hva, hvl
        [1,2,4], ... hba, hva, hvl
        [1,2,5], ... hba, hva, hvf
        [1,3,4], ... hba, hrl, hvl
        [1,3,5], ... hba, hrl, hvf
        ...
        [1,4,5], ... hba, hvl, hvf
        [2,3,4], ... hva, hrl, hvl
        [2,3,5], ... hva, hrl, hvf
        [2,4,5], ... hva, hvl, hvf
        [3,4,5]};  % hrl, hvl, hvf


grps = {[1,3,4],[1,3,5],[1,4,5]};

% labels    'hbang' 'hroll'    'hrvl'    'fhrvf'
grps = {1,2,3}

g = 1;
figure,for x = 1:5,for y = 1:5,subplot2(5,5,y,x);imagesc(mJpdf(:,phzOrder,x,y,1,2,1)');end;end
figure,for x = 1:5,for y = 1:5,subplot2(5,5,y,x);imagesc(mean(mJpdf(:,phzOrder,x,y,3,2,g,:),7)');end;end

figure,imagesc(sq(mean(mmJpdf(2,2,:,:,5,2),2))')

p = phzOrder(2);
e = 2;
g = 4;
figure,
cmap = cool(5);
for z = 1:5
    subplot(1,5,z);
    hold('on');
    for m = 1:3
        for x = 1:5
            plot(sq(mmJpdf(m,p,x,:,z,e,g,:)),'-','Color',cmap(x,:));
        end
    end
    ylim([-300,300]);
end

% STGRPS    [1×4 double]    [4]    [6]    [3]    [5]    [9]    [10]
%%%<<< Display mmJpdf
k = 0;
for s = [1,2,3,4,5];
    k = k+1;
for g = 1:3; 
    for e = 1:2;
        hfig = figure();
        sax = reshape(tight_subplot(8,5,0.005,0.1,0.1),[5,8])';
        for p = 1:8, 
            for z = 1:5, 
                axes(sax(p,z));
                nmask = sq(double(mean(smplCnt(:,z,:,e,g,s,:),7)>500));
                nmask(~nmask) = nan;
                if e==1,
                    ca = [-60,120];
                else
                    ca = [-60, 60];
                end

                imagescnan({vbc{grps{g}(1)},vbc{grps{g}(3)},sq(mean(mmJpdf(2,phzOrder(p),:,z,:,e,g,s,:),9))'.*nmask'},...
                           ca,'linear',false,'colorMap',@jet);
                %imagesc(sq(std(mmJpdf(2,phzOrder(p),:,z,:,e,g,s,:),[],9))');                
                %pcolor(vbc{grps{g}(1)},vbc{grps{g}(3)},sq(mean(mmJpdf(2,phzOrder(p),:,z,:,e,g,s,:),9))');
                axis('xy');
% $$$                 if e==1,
% $$$                     caxis([-20,20]);
% $$$                 else
% $$$                     caxis([-20,20]);
% $$$                 end
                sax(p,z).XTick = [];
                sax(p,z).YTick = [];                        
                if p==1 & z==3,
                    title(vlb{grps{g}(2)});
                end
                if p==8 & z==1,
                    ylabel(vlb{grps{g}(3)});
                    xlabel(vlb{grps{g}(1)});
                end
                
            end            
        end        
        colormap('jet');
        hfig.Units = 'Centimeters';
        hfig.Position = [0+(e-1)*5+10*double(g>4)+10*double(g>8)+(k-1)*2*5,1+(mod(g-1,4))*10+5,5,8];
    end
end
end
%%%>>>

figure();
hold('on');
for x = 
for y = 1:11,    
% $$$     plot([-log10(abs(vbc{5}(vbc{5}<0))),0,log10(vbc{5}(vbc{5}>0))],sq(mean(mmJpdf(2,phzOrder(7),y,6,:,e,g,: ...
% $$$                                                       ),8)),'-+','Color',cmap(y,:));
    plot(vbc{5},sq(mean(mmJpdf(2,phzOrder(7),x,6,:,e,g,:),8)),'-+','Color',cmap(x,:));
end


g = 3;
e = 1;
figure();
cmap = cool(11);
sax = reshape(tight_subplot(8,11,0.005,0.1,0.1),[11,8])';
for p = 1:8
    for y = 1:11,
        axes(sax(p,y));
        %subplot(8,11,y+(p-1)*11);
        hold('on');
        for m = 2,%1:3
            for x = 1:11,
                plot(vbc{grps{g}(3)},sq(mean(mmJpdf(m,phzOrder(p),x,y,:,e,g,:),8)),'-+','Color',cmap(x,:));
% $$$                 plot([-log2(abs(vbc{grps{g}(3)}(vbc{grps{g}(3)}<0))),0,log2(vbc{grps{g}(3)}(vbc{grps{g}(3)}>0))],...
% $$$                      sq(mean(mmJpdf(m,phzOrder(p),x,y,:,e,g,:),8)),'-+','Color',cmap(x,:));
            end
        end
        %xlim([0,12]);
        xlim(vbc{grps{g}(3)}([1,end])+[-5,5]);
        ylim([-30,100]);
        %ylim([-60,60]);        
        Lines([],0,'k');
        Lines(0,[],'k');
    end
end




%%%<<< Compute state residuals 
mmmJpdf = mean(mmJpdf,8);
for g = 1:numel(grps),

    xi = vbi{grps{g}(1)};
    xb = vbc{grps{g}(1)};
    xc = numel(xb);
    yi = vbi{grps{g}(2)};
    yb = vbc{grps{g}(2)};
    yc = numel(yb);
    zi = vbi{grps{g}(3)};
    zb = vbc{grps{g}(3)};
    zc = numel(zb);

    shifts = 0:8:2^8;
    tsmMask = smMask(ind);
    for e = 1:2;
        tferr = lfErrInd{e}(ind);
        tind = ~isnan(tferr+tdphz);
        tind(tind) = logical(tferr(tind)+tdphz(tind));

        % compute JPDF(FERROR, PHZ | HBA, HRVL)
        for x = 1:xc,
            for y = 1:yc,
                for z = 1:zc,
                    disp(['g: ',num2str(g),'   x: ',num2str(x),'  y: ',num2str(y),'  z: ',num2str(z),'  e: ',num2str(e)]);
                    tic
                    for i = 1:100,

                        indb =   (xi==x-1|xi==x|xi==x+1) ...
                               & (yi==y-1|yi==y|yi==y+1) ...
                               & (zi==z-1|zi==z|zi==z+1) ...                                  
                               & circshift(tsmMask, randsample(shifts,1)) ...
                               & tind;
                        indb(randn([sum(indb),1]) > 0) = false;
                        smplCnt(x,y,z,e,g,i) = sum(indb);                        
                        ttferr = tferr(indb);
                        ttphz = tphz(indb);
                        for p = 1:8,    
                             rmJpdf( = prctile((ttferr(ttpbi==p)-mmmJpdf(2,p,x,y,z,e,g))
                        end
                        
                        
                        %mJpdf(:,:,x,y,z,e,g,i) = full(sparse(, tdphz(indb), 1, numFErBins, numPhzBins));
                    end
                    toc
                end
            end
        end
    end
end


%%%>>>







    for vi = 1;
% plot JPDF(FERROR, PHZ | HBA, HRVL)
for e = 1:2;
    hfig = figure();
    hfig.Units = 'centimeters';
    hfig.Position = [(vi-1)*10,16+(e-1)*7,10,numel(hrollBinCenters)];

    sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
    for x = 1:edx; for y = 1:edy;
            axes(sax(edy-y+1,x));               
            hold('on');
            imagesc(ferrorBinEdges{e},                                           ...
                    [phzBins,phzBins+2*pi],                                      ...
                    imgaussfilt(repmat(mirHbaHvlJpdf(:,:,x,y,e),[1,2]),[2,0.1])');
            axis('tight');        axis('xy');
% $$$             text(sax(f,b), -220, -pi,...
% $$$                  ['hba:',num2str(round(hbangBinCenters(b),2)),...
% $$$                   ' hvl:',num2str(round(hrvlBinCenters(f),2))],...
% $$$                  'Color','w',...
% $$$                  'Rotation',90);
    end;end        
    ForAllSubplots('ylim([0,2*pi]);'); 
    if e==1,
        ForAllSubplots('xlim([-200,300]);'); 
    else,
        ForAllSubplots('xlim([-250,250]);');         
    end
    ForAllSubplots('colormap(''jet'');');
    ForAllSubplots('Lines(0,[],''k'');');
    ForAllSubplots('Lines([],0,''k'');');  
    ForAllSubplots('Lines([],pi,''m'');');
    ForAllSubplots('Lines([],2*pi,''k'');');
end

jmask = ones([edx,edy]);
jmask(smplCnt<1e3) = nan;
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [(vi-1)*6,-5,6,numel(hrollBinCenters)*8/3+2];
for p = 1:8,
    subplot2(8,2,p,1);
    imagesc(hbangBinCenters,hrollBinCenters,(sq(mmirHbaHvlJpdf(2,phzOrder(p),:,:,1)).*jmask)');  colormap('jet');  caxis([-60,120]);  axis('xy');
    subplot2(8,2,p,2);
    imagesc(hbangBinCenters,hrollBinCenters,(sq(mmirHbaHvlJpdf(2,phzOrder(p),:,:,2)).*jmask)');  colormap('jet');  caxis([-75,75]);  axis('xy');
end 
end
for j = (25:42);figure(j);end

figure,
hold('on');
cmap = cool(edx);
for a = 1:edx
plot(sq((mean(sq(mmirHbaHvlJpdf(1,phzOrder([2:7]),a,:,2)))))','Color',cmap(a,:))
plot(sq((mean(sq(mmirHbaHvlJpdf(2,phzOrder([2:7]),a,:,2)))))','Color',cmap(a,:))
plot(sq((mean(sq(mmirHbaHvlJpdf(3,phzOrder([2:7]),a,:,2)))))','Color',cmap(a,:))
end

figure();
hold('on');
for a = 1:4
plot(sq((diff(sq(mmirHbaHvlJpdf(2,phzOrder([1,7]),a,:,2)))))','Color',cmap(a,:))
end

figure();
hold('on');
plot(sq((mean(sq(mmirHbaHvlJpdf(1,phzOrder([2:6]),:,:,2)))))'-sq((mean(sq(mmirHbaHvlJpdf(2,phzOrder([2:6]),:,:,2)))))','r')
plot(sq((mean(sq(mmirHbaHvlJpdf(2,phzOrder([2:6]),:,:,2)))))'-sq((mean(sq(mmirHbaHvlJpdf(3,phzOrder([2:6]),:,:,2)))))','b')


figure();
imagesc(sq(mmirHbaHvlJpdf(2,phzOrder([2]),:,:,2))'-sq(mmirHbaHvlJpdf(2,phzOrder([6]),:,:,2))');


figure();
hold('on');
for a = 1:4,
plot(reshape(sq(mmirHbaHvlJpdf(2,phzOrder([5]),a,:,2))',[],1),...
     reshape(sq(mmirHbaHvlJpdf(2,phzOrder([2]),a,:,2))'-sq(mmirHbaHvlJpdf(2,phzOrder([5]),a,:,2))',[],1),...
     '.-');
end


p = [2,5];
mang = -sq(atan2(sq(diff(mmirHbaHvlJpdf(:,phzOrder(p),:,:,1),1,2)),...
              sq(diff(mmirHbaHvlJpdf(:,phzOrder(p),:,:,2),1,2))))'-pi/2;
figure,
subplot(121);
imagescnan({hbangBinCenters,...
           hrvlBinCenters,...
           (mang'.*jmask)'},...
           [-1,1],...
            [],...
            true,...           
           'colorMap',@jet);
axis('xy');
subplot(122);
imagescnan({hbangBinCenters,...
            hrvlBinCenters,...
            (sqrt(sum(sq(diff(mmirHbaHvlJpdf(:,phzOrder(p),:,:,:),1,2)).^2,3)).*jmask)'},...
            [40,100],...
            [],...
            true,...
            'colorMap',@jet);            
axis('xy');

sla = cell([1,2]);
%[sla{:}] = meshgrid([-100,0,100],[-100,0,100]);
[sla{:}] = meshgrid([-150,-75,0,75,150],[-150,-75,0,75,150]);
sla = cf(@(s) reshape(s,[],1), sla);


cmap = cool(7);
figure();
hold('on');
for p = 2,
    quiver(sla{[2,1]},...
           reshape(mmirHbaHvlJpdf(:,phzOrder(p),:,:,2),[],1),...
           reshape(mmirHbaHvlJpdf(:,phzOrder(p),:,:,1),[],1),...
           0,...
           'Color',cmap(p,:));
end
af(@(a) set(a,'ShowArrowHead','off'), findobj(gcf,'Type','quiver'));
set(gca,'YTick',[-100,0,100]);
set(gca,'XTick',[-100,0,100]);
grid(gca,'on');
daspect([1,1,1])
%%%>>>



%%%<<< COMPUTE EPP for up to 2 behavioral variables

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)            ...
        & duincI                                        ...
        & dpostI                                        ...        
        & dhdist<320;

figure()
hist(dhz(ind&~ismember(dtind,[3,4,5])),100);

chzdza = -circ_dist(dhzdza,pi/2);

cblen = dblen+20*double(ismember(dtind,[3,4,5]));
chroll = (dhroll-0.26*double(~ismember(dtind,[3,4,5])));
chbang = -(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5])));
chvang = dhvang;
cfhrvl = dfhrvl;
lfErr = ferr([1:2]);
lfErr = {derrlns,derrlts};
v = cfhrvl;
a = chbang;
chbang  ( a<0 ) = -chbang   ( a<0 );
lfErr{2}( a<0 ) = -lfErr{2} ( a<0 );
cfhrvl  ( a<0 ) = -cfhrvl   ( a<0 );
chroll ( a<0 ) = -chroll  ( a<0 );
chvang  ( a<0 ) = -chvang   ( a<0 );
chzdza  ( a<0 ) = -chzdza   ( a<0 );
% $$$ cfhrvl( a<0 & v>0) = -cfhrvl( a<0 & v>0);
% $$$ cfhrvl( a<0 & v<0) = -cfhrvl( a<0 & v<0);
clear('a');


nbins = 5;
vlb = {};
vdc = {};
vnm = {};
vbe = {};
vbc = {};
vun = {};
%vbi: see below

% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
%vbe{end+1} = [0,0.12,0.26,0.44,0.7,1.2];
switch nbins
  case 3, vbe{end+1} = [0,0.4,0.8,1.2];
  case 5,  vbe{end+1} = linspace(0,1.2,6);
  case 7,  vbe{end+1} = linspace(0,1.2,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';

% 2
vlb{end+1} = 'hva';
vdc{end+1} = 'head anglular velocity';
vnm{end+1} = 'chvang';
switch nbins
  case 3,  vbe{end+1} = [-0.2,-0.012,0.012,0.2];
  case 5,  vbe{end+1} = [-0.2,-0.06,-0.012,0.012,0.06,0.2];    
  case 7,  vbe{end+1} = [-0.2,-0.1,-0.06,-0.012,0.012,0.06,0.1,0.2];
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad/sec';

% 3
vlb{end+1} = 'hrl';
vdc{end+1} = 'head roll';
vnm{end+1} = 'chroll';
switch nbins
  case 3,  vbe{end+1} = [-0.5,-0.1,0.1,0.5];
  case 5,  vbe{end+1} = linspace(-0.5,0.5,6);
  case 7,  vbe{end+1} = linspace(-0.5,0.5,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';

% 4
vlb{end+1} = 'hvl';
vdc{end+1} = 'lateral head speed';
vnm{end+1} = 'cfhrvl';
switch nbins
  case 3,  vbe{end+1} = [-80,-5,5,80];    
  case 5,  vbe{end+1} = [-40,-16,-4,4,16,40];
  case 7,  vbe{end+1} = [-36,-24,-12,-2,2,12,24,36];
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';

% 5
vlb{end+1} = 'hvf';
vdc{end+1} = 'forward head speed';
vnm{end+1} = 'dfhrvf';
switch nbins
  case 3,  vbe{end+1} = [-4,4,12,64];
  case 5,  vbe{end+1} = [-4,4,12,24,36,64];    
  case 7,  vbe{end+1} = [-4,0,4,10,20,30,40,80];
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';

% 6
vlb{end+1} = 'hbd';
vdc{end+1} = 'head-body distance';
vnm{end+1} = 'cblen';
switch nbins
  case 3,  vbe{end+1} = [100,120,140,160];
  case 5,  vbe{end+1} = linspace(100,160,6);
  case 7,  vbe{end+1} = linspace(100,170,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'mm';

% 7
vlb{end+1} = 'hz';
vdc{end+1} = 'head height';
vnm{end+1} = 'dhz';
switch nbins
  case 3,  vbe{end+1} = [40,60,80,100];
  case 5,  vbe{end+1} = linspace(40,100,6);
  case 7,  vbe{end+1} = linspace(40,100,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm';

% 8
vlb{end+1} = 'hbp';
vdc{end+1} = 'head-body pitch';
vnm{end+1} = 'dfet';
switch nbins
  case 3,  vbe{end+1} = [-1.4,-0.8,-0.2,0.2];
  case 5,  vbe{end+1} = linspace(-1.4,0.2,6);
  case 7,  vbe{end+1} = linspace(-1.4,0.2,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';

% 9
vlb{end+1} = 'hvd';
vdc{end+1} = 'head-body distance velocity';
vnm{end+1} = 'dhvd';
switch nbins
  case 3,  vbe{end+1} = [-80,-5,5,80];
  case 5,  vbe{end+1} = linspace(-80,80,6);
  case 7,  vbe{end+1} = linspace(-80,80,8);
end
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';


vlb{9} = 'hda';
vdc{9} = 'hrz-drz angle';
vnm{9} = 'chzdza';
switch nbins
  case 3,  vbe{9} = [-pi,-pi/3,pi/3,pi];
  case 5,  vbe{9} = linspace(-pi,pi,6);
  case 7,  vbe{9} = linspace(-pi,pi,8);
end
vbc{9} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{9} = 'rad';






stgrps = {[9:11],4,6,3,5,9,10,11};
stlbls = {'all','HP','LP','HW','LW','T','W','P'};
nvars = numel(vlb);

tag = DataHash({vlb,stgrps,nbins,nvars,stlbls});
filepath = fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_2d_',tag,'.mat']);

if exist(filepath,'file'),
    load(filepath);
else
    vbi = {};
    sJpdf =  zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),nbins, nbins, 2,nvars,nvars,numel(stgrps)]);
    smJpdf = zeros([3,numel(phzBinCenters),nbins, nbins, 2,nvars,nvars,numel(stgrps)]);
    sCnt =  zeros([nbins, nbins, 2,nvars,nvars,numel(stgrps)]);
    for s = 1:numel(stgrps),
        ind = logical(dstcm(:,1))                             ...
              & any(logical(dstcm(:,stgrps{s})),2)            ...
              & duincI                                        ...
              & dpostI                                        ...        
              & dhdist<320;% & chbang <0.25;
        for v = 1:numel(vbe),
            vbi{v} = discretize(eval([vnm{v},'(ind)']),vbe{v});
        end
        tdphz = dphz(ind);
        for g = 1:nvars,
            for h = g+1:nvars,
                xi = vbi{g};        xb = vbc{g};        xc = numel(xb);
                yi = vbi{h};        yb = vbc{h};        yc = numel(yb);
                for e = 1:2;
                    tferr = lfErr{e}(ind);
                    % compute JPDF(FERROR, PHZ | HBA, HRVL)
                    for x = 1:xc,
                        for y = 1:yc,
+                            disp(['s: ',num2str(s),'   g: ',num2str(g),'   h: ',num2str(h),'   x: ',num2str(x),'  y: ',num2str(y),'  e: ',num2str(e)]);
                            indb = xi==x & yi==y;
                            sCnt(x,y,e,g,h,s) = sum(indb);
                            sJpdf(:,:,x,y,e,g,h,s) = hist2([tferr(indb),tdphz(indb)],...
                                                           ferrorBinEdges{e},...
                                                           phzBins);
                        end
                    end
                end
                
            end
        end
    end

    for s = 1:numel(stgrps),
        for g = 1:nvars,
            for h = g+1:nvars,        
                for x = 1:xc,
                    tic
                    for y = 1:yc
                        for p = 1:numel(phzBinCenters),
                            for e = 1:2,
                                try
                                    [c, index] = unique(cumsum(sJpdf(:,p,x,y,e,g,h,s)./sum(sJpdf(:,p,x,y,e,g,h,s))),'last');
                                    index(c==0|c==1) = [];
                                    c(c==0|c==1) = [];
                                    smJpdf(:,p,x,y,e,g,h,s) = interp1(c,ferrorBinCenters{e}(index),[0.25,0.5,0.75]);
                                end
                            end
                        end
                    end
                    toc
                end
            end
        end
    end

    save(fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_2d_',tag,'.mat']),...
         'smJpdf',                                                                           ...
         'sJpdf',                                                                            ...
         'sCnt',                                                                             ...
         'vlb',                                                                              ...
         'vdc',                                                                              ...
         'vbe',                                                                              ...
         'vbc',                                                                              ...
         'vnm',                                                                              ...     
         'vun',                                                                              ...          
         'phzBins',                                                                          ...
         'stgrps',                                                                           ...
         '-v7.3'                                                                             ...
         );
end

%%%>>>

pbc = phzBinCenters+2*pi*double(phzBinCenters<0);


%1 hba
%2 hva
%3 hrl
%4 hvl
%5 hvf
%6 hbd
%7 hz
%8 hbp



s = 6;
e = 2;
g =4;
h = 5;
figure();
sax = reshape(tight_subplot(5,5,[0.01,0.01],0.1,0.1)',[5,5]);
for x = 1:size(sJpdf,3);
    for y = 1:size(sJpdf,4);
        axes(sax(x,6-y));
        imagesc(ferrorBinCenters{1},pbc(phzOrder),imgaussfilt(sq(sJpdf(:,phzOrder,x,y,e,g,h,s))',[0.1,3]));
        %imagesc(ferrorBinCenters{1},pbc(phzOrder),sq(sJpdf(:,phzOrder,x,y,e,g,h,s))');        
        axis('xy');
        if e==1,
            xlim([-200,300]);
        else
            xlim([-300,300]);
        end
        Lines(0,[],'k');
        Lines(100,[],'k');
        Lines(-100,[],'k');
        sax(x,6-y).XTick = [];
        sax(x,6-y).YTick = [];      
        if x ==1 && y == 3,
            ylabel(vlb{h});
        end
        if x == 3 && y == 1,
            xlabel(vlb{g});
        end
    end
end
colormap('jet');



s = 1;
ind = logical(dstcm(:,1))                             ...
      & any(logical(dstcm(:,stgrps{s})),2)            ...
      & duincI                                        ...
      & dpostI                                        ...        
      & dhdist<380;% & chbang <0.25;
for g = 1:nvars,
    for h = 1:nvars,
        subplot2(nvars,nvars,h,g);
        hist2([eval([vnm{g},'(ind)']),eval([vnm{h},'(ind)'])],...
              linspace([vbe{g}([1,end]),50]),...
              linspace([vbe{h}([1,end]),50]));
        Lines(vbe{g},[],'k');
        Lines([],vbe{h},'k');
        title(stlbls{s});
    end
end

%%%<<< DISPLAY the jdpf of space and median EPP F and L
g = 1;
h = 4;
figure();
k = 0;
for s = [1:8];
    %for s = [1,6,7,8];
    k = k +1;
    subplot(1,8,k);
    ind = logical(dstcm(:,1))                             ...
          & any(logical(dstcm(:,stgrps{s})),2)            ...
          & duincI                                        ...
          & dpostI                                        ...        
          & dhdist<380;% & chbang <0.25;
    hist2([eval([vnm{g},'(ind)']),eval([vnm{h},'(ind)'])],...
          linspace([vbe{g}([1,end]),50]),...
          linspace([vbe{h}([1,end]),50]));
    Lines(vbe{g},[],'k');
    Lines([],vbe{h},'k');
    title(stlbls{s});
end



k = 0;
for s = [7];
    %for s = [1,6,7,8];
    k = k+1;
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [(k-1).*6,5,6,24];
for p = 1:8,
    for e = 1:2
        subplot2(8,2,p,e);
        nmask = sq(double(sCnt(:,:,2,g,h,s)>10000));
        nmask(~nmask) = nan;
        if e==1,
            ca = [-60,120];
        else
            ca = [-60, 60];
        end
% $$$         imagescnan({vbc{grps{g}(1)},vbc{grps{g}(2)},sq(smJpdf(2,phzOrder(p),:,:,e,g,s,:))'.*nmask'},...
% $$$                    ca,'linear',false,'colorMap',@jet);
% $$$         axis('xy');
        imagesc(vbc{g},vbc{h},sq(smJpdf(2,phzOrder(p),:,:,e,g,h,s,:))'.*nmask');
        colormap('jet');
        caxis(ca)
        axis('xy');
    
        if e==1&&p==1,
            title(stlbls{s});
        end
    end
end
end

%%%>>>


p = 2;
figure,imagesc(sqrt(sum(sq(diff(smJpdf(2,phzOrder([2,6]),:,:,:,g,h,s,:),1,2)).^2,3))'.*nmask');
axis('xy');
colormap('jet');
caxis([0,100])



ind = logical(dstcm(:,1))                             ...
      & duincI                                        ...
      & dpostI                                        ...        
      & dhdist<320;% & chbang <0.25;


tp = ThreshCross(circshift(dstcm(:,11),64)&dstcm(:,9),0.5,32);

ssegs = GetSegs(dstcm(:,[9,11]),tp(:,1)-4000,8000);
tsegs = GetSegs(lfErr{2},tp(:,1)-4000,8000);
vsegs = GetSegs(cfhrvl,tp(:,1)-4000,8000);
tinds = GetSegs(ind,tp(:,1)-4000,8000,false);
asegs = GetSegs(chbang,tp(:,1)-4000,8000);

tsegs(~tinds)=nan;
ssegs(:,sum(~isnan(tsegs))<4000,:) = [];
asegs(:,sum(~isnan(tsegs))<4000) = [];
vsegs(:,sum(~isnan(tsegs))<4000) = [];
tsegs(:,sum(~isnan(tsegs))<4000) = [];

figure,imagesc(tsegs(3:8:end,:)')


pbt = mode(sum(ssegs(2000:4000,:,:),3))==11;
%pbt = true;
tf = 4240;
tt = 4000;
%tt = 4480;
%tt = 3760;
aind = {};
aind{1} = asegs(tt,:) > 0.0 & asegs(tt,:) < 0.4;% & pbt & mean(abs(vsegs(tt-100:tt+100,:)))>10;
aind{2} = asegs(tf,:) > 0.4 & asegs(tf,:) < 0.8;% & pbt & mean(abs(vsegs(tt-100:tt+100,:)))>10;
aind{3} = asegs(tt,:) > 0.8 & asegs(tt,:) < 1.4;% & pbt & vsegs(tf,:)<0;
cmap = cool(6);
figure();
for a = 1:3,
    subplot(2,3,a);hold('on');
    for c = 1:6,
        plot(linspace(-2,2,1000),median(tsegs(phzOrder(c):8:end,aind{a}),2,'omitnan'),'Color',cmap(7-c,:));
    end
    Lines(0,[],'k');
    Lines([-0.15,0.15],[],'r');    
    
    subplot(2,3,a+3);hold('on');
    for c = 1:6,
        plot(linspace(-2,2,1000),std(tsegs(phzOrder(c):8:end,aind{a}),[],2,'omitnan'),'Color',cmap(7-c,:));
    end
    
    Lines(0,[],'k');
    Lines([-0.15,0.15],[],'r');    
end

figure,plot(mean(sqrt(sum(sq(diff(smJpdf(2,phzOrder([2,5]),:,:,:,g,h,s,:),1,2)).^2,3)).*nmask));

%%%<<< REPORT FIGURE: jdpf of space and median EPP F and L
hfig = figure();
hfig.Units = 'Centimeters';
hfig.Position = [5,5,29.7,21];
pause(0.1);
pause(0.1);

g = 1;
h = 6;

clf();
k = 0;
try,
    delete(sp);
    delete(tp);
end
sp = gobjects([0,1]);
tp = gobjects([0,1]);
cax = gobjects([0,1]);
for s = [1:8];
    %for s = [1,6,7,8];
    k = k +1;
    tp(end+1) = subplot2(11,16,[1,2],[1:2]+(k-1).*2);
    ind = logical(dstcm(:,1))                             ...
          & any(logical(dstcm(:,stgrps{s})),2)            ...
          & duincI                                        ...
          & dpostI                                        ...        
          & dhdist<380;% & chbang <0.25;
    hist2([eval([vnm{g},'(ind)']),eval([vnm{h},'(ind)'])],...
          linspace([vbe{g}([1,end]),50]),...
          linspace([vbe{h}([1,end]),50]));
    Lines(vbe{g},[],'k');
    Lines([],vbe{h},'k');
    title(stlbls{s});
    if k ~= 1;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
    else
        ylabel(vun{h});
        xlabel(vun{g});
        title([stlbls{s},' - ',vlb{g},' VS ',vlb{h}]);
    end
    
    for p = 1:8,
        for e = 1:2
            sp(end+1) = subplot2(11,16,p+3,2*(k-1)+e);
            nmask = sq(double(sCnt(:,:,2,g,h,s)>10000));
            nmask(~nmask) = nan;
            if e==1,
                ca = [-60,120];
            else
                ca = [-60, 60];
            end
            [~,cax(end+1) ] = imagescnan({vbc{g},vbc{h},sq(smJpdf(2,phzOrder(p),:,:,e,g,h,s))'.*nmask'},...
                                  ca,'linear',k==8&p==8,'colorMap',@jet);
            axis('xy');
            if e==1&&p==1,
                title('FWD');
            elseif e==2&&p==1,
                title('LAT');
            end
            if k>1||e==2,
                sp(end).YTick = [];
            else
                ylabel(vun{h});                
            end
            
            if p < 8 || e==2
                sp(end).XTick = [];
            else
                xlabel(vun{g});
                sp(end).XTick = vbc{g}([1,end]);
            end
            
            if k == 8 && e==1
                xlabel('');
            end
        end
    end
end


af(@(s) set(s,'Position',get(s,'Position').*[1,1,0,0]+get(sp(1),'Position').*[0,0,1,1]), sp(end-1:end));

af(@(s) set(s,'Units','Centimeters'), sp);
af(@(s) set(s,'Position',get(s,'Position')+[-0.25,-0.5,0.25,0.5]), sp);
af(@(s) set(s,'Units','Centimeters'), tp);
af(@(s) set(s,'Position',get(s,'Position')+[-0.25,-1,0,0]), tp);

af(@(c) set(c,'Units','Centimeters'),                  cax(end-1:end));
af(@(c) set(c,'Position',get(c,'Position')-[0,1,0,0]), cax(end-1:end));
af(@(c,s) set(c,'Position',s.Position.*[1,1,1,0]+[0,-1.5,0,0.25]), cax(end-1:end),sp(end-1:end));

af(@(c) view(c,[90,90]), cax(end-1:end));
cax(end-1).XTickLabels = [-50,100];
cax(end).XTickLabels = [-60,60];
cax(end-1).Position(1) = cax(end-1).Position(1) - 0.25;
cax(end).Position(1) = cax(end).Position(1) + 0.25;
sp(end-1).XTick = [];
ylabel(cax(end-1),'FWD')
ylabel(cax(end),'LAT')

fax = axes('Position',[0,0,1,1],'Visible','off');
text(0.1,0.95,{[vlb{g},': ',vdc{g}],[vlb{h},': ',vdc{h}]});

filePath = fullfile(MTA_PROJECT_REPORT_PATH,'egocentricPP',['MjgER2016_hbaCorrected_2d_',tag,'_',vlb{g},'_',vlb{h}]);
print(hfig,'-dpng', [filePath,'.png']);
%print(hfig,'-depsc',[filePath,'.eps']);



s = 7;
g = 5;
h = 6;
e = 1;
n = 3;
d = 5;
cmap = cool(nbins);
hfig = figure();
hold('on');
for c = 2:nbins,
    plot(sq(mean(smJpdf(2,phzOrder([4:5]),c,n:d,e,g,h,s))),...
         sq(diff(smJpdf(2,phzOrder([6,2]),c,n:d,e,g,h,s))),...
         '-+','Color',cmap(c,:));
end
legend({'-4 to 4 cm/s','4 to 12 cm/s','12 to 24 cm/s','24 to 36 cm/s','36 to 64 cm/s'});
title({'Theta trough FEPP vs theta cycle FEPP span','WALK; HBD=[130*,142,154]'});
for c = 2:nbins,
    plot(sq(mean(smJpdf(2,phzOrder([4:5]),c,n,e,g,h,s))),...
         sq(diff(smJpdf(2,phzOrder([6,2]),c,n,e,g,h,s))),...
         '*','Color',cmap(c,:));
end
ylabel('FEPP span (mm)');
xlabel('Theta Trough FEPP (mm)');
filePath = fullfile(MTA_PROJECT_REPORT_PATH,'egocentricPP',['hvf_hbd_walk_FEPP_Theta_Trough_',tag]);
print(hfig,'-dpng', [filePath,'.png']);


clf();
hold('on');
for c = 1:nbins,
    plot(sq(mean(smJpdf(2,phzOrder([2:3]),c,n:d,e,g,h,s))),...
         sq(diff(smJpdf(2,phzOrder([6,2]),c,n:d,e,g,h,s))),...
         '-+','Color',cmap(c,:));
end
lax = legend({'-4 to 4 cm/s','4 to 12 cm/s','12 to 24 cm/s','24 to 36 cm/s','36 to 64 cm/s'},'Location','southeast');
title({'Theta Ascent FEPP vs theta cycle FEPP span','WALK; HBD=[130*,142,154]'});
for c = 1:nbins,
    plot(sq(mean(smJpdf(2,phzOrder([2:3]),c,n,e,g,h,s))),...
         sq(diff(smJpdf(2,phzOrder([6,2]),c,n,e,g,h,s))),...
         '*','Color',cmap(c,:));
end
ylabel('FEPP span (mm)');
xlabel('Theta Ascent FEPP (mm)');
filePath = fullfile(MTA_PROJECT_REPORT_PATH,'egocentricPP',['hvf_hbd_walk_FEPP_Theta_Ascent',tag]);
print(hfig,'-dpng', [filePath,'.png']);


clf();
hold('on');
for c = 1:nbins,
    plot(sq(mean(smJpdf(2,phzOrder([6:7]),c,n:d,e,g,h,s))),...
         sq(diff(smJpdf(2,phzOrder([6,2]),c,n:d,e,g,h,s))),...
         '-+','Color',cmap(c,:));
end
lax = legend({'-4 to 4 cm/s','4 to 12 cm/s','12 to 24 cm/s','24 to 36 cm/s','36 to 64 cm/s'},'Location','southeast');
title({'Theta Descent FEPP vs theta cycle FEPP span','WALK; HBD=[130*,142,154]'});
for c = 1:nbins,
    plot(sq(mean(smJpdf(2,phzOrder([6:7]),c,n,e,g,h,s))),...
         sq(diff(smJpdf(2,phzOrder([6,2]),c,n,e,g,h,s))),...
         '*','Color',cmap(c,:));
end
ylabel('FEPP span (mm)');
xlabel('Theta Descent FEPP (mm)');
filePath = fullfile(MTA_PROJECT_REPORT_PATH,'egocentricPP',['hvf_hbd_walk_FEPP_Theta_Descent',tag]);
print(hfig,'-dpng', [filePath,'.png']);




%%%>>>


%%%<<< Indepth 2d 

vlb = {};
vdc = {};
vnm = {};
vbe = {};
vbc = {};
vun = {};
%vbi: see below


%%%<<< HBA x HVF
% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
vbe{end+1} = [0,0.4,0.8,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';
%2
vlb{end+1} = 'hvf';
vdc{end+1} = 'forward head speed';
vnm{end+1} = 'dfhrvf';
%vbe{end+1} = [-10,10,30,50,80];
%vbe{end+1} = [-4,4,12,18,24,30,36,42,48,54,60,66,72,78];
%vbe{end+1} = [-4,4,10,20,30,40,50,60,70,80];
vbe{end+1} = [-4,4,10,25,40,55,70,85];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';
%%%>>>

%%%<<< HBA x HVF
% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
vbe{end+1} = [0,0.4,0.8,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';
%2
vlb{end+1} = 'hvl';
vdc{end+1} = 'lateral head speed';
vnm{end+1} = 'dfhrvl';
%vbe{end+1} = [-10,10,30,50,80];
%vbe{end+1} = [-4,4,12,18,24,30,36,42,48,54,60,66,72,78];
%vbe{end+1} = [-4,4,10,20,30,40,50,60,70,80];
vbe{end+1} = [-60,-40,-25,-10,-4,4,10,25,40,60];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'cm/s';
%%%>>>


%%%<<< HBA x HBD
% 1
vlb{end+1} = 'hba';
vdc{end+1} = 'head-body angle';
vnm{end+1} = 'chbang';
vbe{end+1} = [0,0.4,0.8,1.2];
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'rad';
% 2
vlb{end+1} = 'hbd';
vdc{end+1} = 'head-body distance';
vnm{end+1} = 'cblen';
vbe{end+1} = linspace(110,160,9);
vbc{end+1} = mean([vbe{end}(2:end); vbe{end}(1:end-1)]);
vun{end+1} = 'mm';
%%%>>>



stgrps = {9:11,2,9,10,11};
stlbls = {'ALL','Rear','Turn','Walk','Pause'};
nIter = 100;
shifts = 0:8:2^8;
tag = DataHash({vlb,stgrps,nbins,nvars,stlbls});
filepath = fullfile(MTA_PROJECT_PATH,'analysis',['MjgER2016_req20191104_jpdf_hbaCorrected_2d_ID_',tag,'.mat']);

if exist(filepath,'file'),
    load(filepath);
else

    vbi = {};
    sJpdf =  zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2, numel(stgrps),nIter]);
    smJpdf = zeros([3,numel(phzBinCenters),numel(vbc{1}),numel(vbc{2}), 2,numel(stgrps),nIter]);
    sCnt =  zeros([numel(vbc{1}),numel(vbc{2}), 2,numel(stgrps),nIter]);
    for s = 4%:numel(stgrps),
        ind = logical(dstcm(:,1))                             ...
              & any(logical(dstcm(:,stgrps{s})),2)            ...
              & duincI                                        ...
              & dpostI                                        ...        
              & dhdist<320;% & chbang<0.4;
        for v = 1:numel(vbe),
            vbi{v} = discretize(eval([vnm{v},'(ind)']),vbe{v});
        end
        tdphz = dphz(ind);
        tsmMask = smMask(ind);
        for g = 1,
            for h = 2,
                xi = vbi{g};        xb = vbc{g};        xc = numel(xb);
                yi = vbi{h};        yb = vbc{h};        yc = numel(yb);
                for e = 1:2,
                    tferr = lfErr{e}(ind);
                    % compute JPDF(FERROR, PHZ | HBA, HRVL)
                    for x = 1:xc,
                        for y = 1:yc,
                            disp(['s: ',num2str(s),'   x: ',num2str(x),'  y: ',num2str(y),'  e: ',num2str(e)]);
                            for i = 1:nIter,
                                indb = xi==x ...
                                       & yi==y ...
                                       & circshift(tsmMask, randsample(shifts,1));
                                indb(randn([sum(indb),1]) > 0) = false;
                                sCnt(x,y,e,s,i) = sum(indb);
                                sJpdf(:,:,x,y,e,s,i) = hist2([tferr(indb),tdphz(indb)],...
                                                             ferrorBinEdges{e},...
                                                             phzBins);
                            end
                        end
                    end
                end
                
            end
        end
    end

    for s = 4%:numel(stgrps)
        for x = 1:xc,
            tic
            for y = 1:yc
                for p = 1:numel(phzBinCenters),
                    for e = 1:2,
                        for i = 1:nIter
                        try
                            [c, index] = unique(cumsum(sJpdf(:,p,x,y,e,s,i)./sum(sJpdf(:,p,x,y,e,s,i))),'last');
                            index(c==0|c==1) = [];
                            c(c==0|c==1) = [];
                            smJpdf(:,p,x,y,e,s,i) = interp1(c,ferrorBinCenters{e}(index),[0.25,0.5,0.75]);
                        end
                        end
                    end
                end
            end
            toc
        end
    end

    save(filepath,...
         'smJpdf',                                                                           ...
         'sJpdf',                                                                            ...
         'sCnt',                                                                             ...
         'vlb',                                                                              ...
         'vdc',                                                                              ...
         'vbe',                                                                              ...
         'vbc',                                                                              ...
         'vnm',                                                                              ...     
         'vun',                                                                              ...          
         'phzBins',                                                                          ...
         'stgrps',                                                                           ...
         '-v7.3'                                                                             ...
         );
end




k = 0;
%for s = [1:5];
for s = [4];
    %for s = [1,6,7,8];
    k = k+1;
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [(k-1).*6,5,6,24];
for p = 1:8,
    for e = 1:2
        subplot2(8,2,p,e);
        nmask = double(sq(mean(sCnt(:,:,e,s,1:nIter),ndims(sCnt)))>300);
        nmask(~nmask) = nan;
        %nmask = 1;
        if e==1,
            ca = [-120,200];
        else
            ca = [-60, 60];
        end
% $$$         imagescnan({vbc{grps{g}(1)},vbc{grps{g}(2)},sq(smJpdf(2,phzOrder(p),:,:,e,g,s,:))'.*nmask'},...
% $$$                    ca,'linear',false,'colorMap',@jet);**
% $$$         axis('xy');
        imagesc(vbc{1},vbc{2},sq(median(smJpdf(2,phzOrder(p),:,:,e,s,:),ndims(smJpdf)))'.*nmask');
        colormap('jet');
        caxis(ca)
        axis('xy');
    
        if e==1&&p==1,
            title(stlbls{s});
        end
    end
end
end

s = 4;
x = 1;
figure();
for y = 1:numel(vbc{2});
    subplot2(numel(vbc{2}),2,y,1);
    imagesc(ferrorBinCenters{1},...
            phzBinCenters(phzOrder)+2*pi*double(phzBinCenters(phzOrder)<0), ...
            imgaussfilt(sq(mean(sJpdf(:,phzOrder,x,y,1,s,:),ndims(sJpdf))),[2,0.1])');
    axis('xy');
    Lines(0,[],'k');
    Lines([],pi,'k');
    xlim([-250,350]);
    subplot2(numel(vbc{2}),2,y,2);
    imagesc(ferrorBinCenters{1},...
            phzBinCenters(phzOrder)+2*pi*double(phzBinCenters(phzOrder)<0), ...
            imgaussfilt(sq(mean(sJpdf(:,phzOrder,x,y,2,s,:),ndims(sJpdf))),[2,0.1])');
    axis('xy');
    Lines(0,[],'k');
    Lines([],pi,'k');
    xlim([-250,250]);
end
colormap('jet');

e = 1;
x = 2;
s = 3;
switch s,
  case 1,
    span = 1:numel(vbc{2});
  case 2,
    span = 2:numel(vbc{2})-3;
  case 3,
    span = 3:numel(vbc{2});
  case 4,
    span = 1:numel(vbc{2})-1;
end

figure();
hold('on');
cmap = cool(numel(vbc{2}));
for y = span,
    for i = 1:nIter,
    plot(sq(smJpdf(2,phzOrder(1:end),x,y,e,s,i)),'-+','Color',cmap(y,:));
    end
end
for y = span,
    plot(sq(mean(smJpdf(2,phzOrder(1:end),x,y,e,s,:),ndims(smJpdf))),'-+','Color','k','LineWidth',2);
end
ylim([-25,135]);








g = 1;
h = 2;
figure();
k = 0;
for s = [1:5]
    %for s = [1,6,7,8];
    k = k +1;
    subplot(1,4,k);
    ind = logical(dstcm(:,1))                             ...
          & any(logical(dstcm(:,stgrps{s})),2)            ...
          & duincI                                        ...
          & dpostI                                        ...        
          & dhdist<380;% & chbang <0.25;
    hist2([eval([vnm{g},'(ind)']),eval([vnm{h},'(ind)'])],...
          linspace([vbe{g}([1,end]),50]),...
          linspace([vbe{h}([1,end]),50]));
    Lines(vbe{g},[],'k');
    Lines([],vbe{h},'k');
    title(stlbls{s});
end

%%%>>>



% BS mJpdf
mJpdf =  zeros([numel(ferrorBinCenters{1}),numel(phzBinCenters),nbins, nbins, 2,numel(grps),numel(stgrps),100]);
smplCnt =  zeros([nbins, nbins, 2,numel(grps),numel(stgrps),100]);
mmJpdf = zeros([3,numel(phzBinCenters),nbins, nbins, 2,numel(grps),numel(stgrps),100]);
for s = 1:numel(stgrps),
    ind = logical(dstcm(:,1))                             ...
          & any(logical(dstcm(:,stgrps{s})),2)            ...
          & duincI                                        ...
          & dpostI                                        ...        
          & dhdist<380;
    vbi{1} = discretize(chbang(ind),vbe{1});
    vbi{2} = discretize(chvang(ind),vbe{2});
    vbi{3} = discretize(chroll(ind),vbe{3});
    vbi{4} = discretize(cfhrvl(ind),vbe{4});
    vbi{5} = discretize(dfhrvf(ind),vbe{5});
    vbi{6} = discretize(cblen(ind),vbe{6});    

    tdphz = dphz(ind);
    
    for g = 1:numel(grps),

        xi = vbi{grps{g}(1)};
        xb = vbc{grps{g}(1)};
        xc = numel(xb);
        yi = vbi{grps{g}(2)};
        yb = vbc{grps{g}(2)};
        yc = numel(yb);

        shifts = 0:8:2^8;
        tsmMask = smMask(ind);
        for e = 1:2;
            tferr = lfErr{e}(ind);
            % compute JPDF(FERROR, PHZ | HBA, HRVL)
            for x = 1:xc,
                for y = 1:yc,
                    disp(['s: ',num2str(s),'   g: ',num2str(g),'   x: ',num2str(x),'  y: ',num2str(y),'  e: ',num2str(e)]);
                    tic
                    for i = 1:100,
                        indb = xi==x & yi==y  &  circshift(tsmMask, randsample(shifts,1));
                        indb(indb) = randn([sum(indb),1]) > 0;
                        smplCnt(x,y,e,g,s,i) = sum(indb);
                        mJpdf(:,:,x,y,e,g,s,i) = hist2([tferr(indb),tdphz(indb)],...
                                                       ferrorBinEdges{e},...
                                                       phzBins);
                    end
                    toc
                end
            end
        end
    end
end


for s = 1:numel(stgrps),
for g = 1:numel(grps)
    for x = 1:xc,
        tic
        for y = 1:yc
            for p = 1:numel(phzBinCenters),
                for e = 1:2,
                    for i = 1:100,
                        try
                            [c, index] = unique(cumsum(mJpdf(:,p,x,y,e,g,s,i)./sum(mJpdf(:,p,x,y,e,g,s,i))),'last');
                            index(c==0|c==1) = [];
                            c(c==0|c==1) = [];
                            mmJpdf(:,p,x,y,e,g,s,i) = interp1(c,ferrorBinCenters{e}(index),[0.25,0.5,0.75]);
                        end
                    end
                end
            end
        end
        toc
    end
end
end

save(fullfile(MTA_PROJECT_PATH,'analysis','MjgER2016_req20191104_jpdf_hbaCorrected_2d.mat'),...
     'mmJpdf',                                                                           ...
     'mJpdf',                                                                            ...
     'vlb',                                                                              ...
     'vbe',                                                                              ...
     'vbc',                                                                              ...
     'phzBins',                                                                          ...
     'grps',                                                                             ...
     'stgrps',                                                                           ...
     '-v7.3'                                                                             ...
);
load(fullfile(MTA_PROJECT_PATH,'analysis','MjgER2016_req20191104_jpdf_hbaCorrected_2d.mat'));



cblen = dblen+20*double(ismember(dtind,[3,4,5]));

s = 7
    ind = logical(dstcm(:,1))                             ...
          & any(logical(dstcm(:,[3,5,4,6])),2)            ...
          & duincI                                        ...
          & dpostI                                        ...        
          & dhdist<380;

    ind = logical(dstcm(:,1))                             ...
          & any(logical(dstcm(:,[5])),2)            ...
          & duincI                                        ...
          & dpostI                                        ...        
          & dhdist<380                                    ...        
          & chbang<0.25;
    
figure();
for p = 1:8,
    subplot2(8,2,p,1);
    idx = ind&dphz==phzBinCenters(phzOrder(p))&smMask;%&dfhrvf>5;
    hist2([ferr{1}(idx),...
           cblen(idx)],...
           linspace(-200,300,30),...
           linspace(110,160,30));
    subplot2(8,2,p,2);    
    hist2([ferr{1}(idx),...
           sqrt(dfhrvf(idx))],...
          linspace(-200,300,30),...
          linspace(1,8,30));
end
colormap('jet');

idx = idx&dfhrvf>1;
figure,plot3(ferr{1}(idx),...
           sqrt(dfhrvf(idx)),...
                        cblen(idx),'.');


figure();
for p = 1:8,
    subplot(8,1,p);
    hist2([ferr{1}(ind&dphz==phzBinCenters(phzOrder(p))),...
           sqrt(dfhrvf(ind&dphz==phzBinCenters(phzOrder(p))))],...
           linspace(-200,300,30),...
           linspace(1,8,30));
end
colormap('jet');

figure();
hist2([log10(dxyvel(ind&dphz==phzBinCenters(3))),dblen(ind&dphz==phzBinCenters(3))],100,100)
figure();
hist2([ferr{1}(ind&dphz==phzBinCenters(3)),log10(dxyvel(ind&dphz==phzBinCenters(3)))],100,100)

figure();
hist2([ferr{1}(ind&dphz==phzBinCenters(3)),dfet(ind&dphz==phzBinCenters(3))],100,100)

figure,hist2([chbang(ind),cblen(ind)],linspace(0,1.5,50),linspace(110,160,50));

figure,hist2([dfet(ind),cblen(ind)],linspace(-1.5,0.25,50),linspace(110,160,50));

p = 2;
idx = ind&dphz==phzBinCenters(phzOrder(p))&smMask;
[R,P] = corrcoef(ferr{1}(idx),dblen(idx))
mdl = glmfit(ferr{1}(idx),dblen(idx))

[B,DEV,STATS] = glmfit( [dfhrvf(idx),cblen(idx)],[dErrlon(idx),);

% g
% 1 [1,4] hba,hvl
% 2 [1,5] hba,hvf
% 3 [1,3] hba,hrl
% 4 [5,6] hvf,bln
% 5 [4,6] hvl,bln
% 6 [1,6] hba,bln

g = 4;
for s = 1:7;
hfig = figure();
hfig.Units = 'centimeters';
hfig.Position = [(s-1).*4,10,4,16];
for p = 1:8,
    for e = 1:2
        subplot2(8,2,p,e);
        nmask = sq(double(median(smplCnt(:,:,2,g,s,:),ndims(smplCnt),'omitnan')>500));
        nmask(~nmask) = nan;
        if e==1,
            ca = [-60,120];
        else
            ca = [-60, 60];
        end
        imagescnan({vbc{grps{g}(1)},vbc{grps{g}(2)},sq(mean(mmJpdf(2,phzOrder(p),:,:,e,g,s,:),ndims(mmJpdf),'omitnan'))'.*nmask'},...
                   ca,'linear',false,'colorMap',@jet);
        axis('xy');
    end
end
end






figure,
for y = 1:7,
    subplot(7,1,8-y);
    imagesc(sq(mean(mJpdf(:,phzOrder,6,y,2,g,s,:),ndims(mJpdf),'omitnan'))');
end

g = 2;
figure();
xc = numel(vbc{grps{g}(1)});
yb = vbc{grps{g}(2)};
cmap = cool(xc);
scl = 10;
sax = reshape(tight_subplot(8,2,[0.01,0.1],0.1,0.1)',[2,8]);
for e = 1:2,
for p = 1:8,    
axes(sax(e,p));
hold('on');
for x = 1:xc
errorbar(yb,sq(mean(mmJpdf(2,phzOrder(p),x,:,e,g,:),7)),sq(std(mmJpdf(2,phzOrder(p),x,:,e,g,:),[],7)./scl),'Color',cmap(x,:));
end
if e==1,
    ylim([-20,90]);
else
    ylim([-20,60]);
end
Lines([],0,'k');
end
end


figure,
sax = reshape(tight_subplot(8,1,[0.01,0.1],0.1,0.1)',[1,8]);
for p = 1:8,    
axes(sax(p));
hold('on');
for x = 1:xc
plot(sq(mean(mmJpdf(2,phzOrder(p),x,:,2,g,:),7)),sq(mean(mmJpdf(2,phzOrder(p),x,:,1,g,:),7)),'Color',cmap(x,:));
end
ylim([-20,90]);
xlim([-20,55]);
Lines([],0,'k');
Lines(0,[],'k');
end


figure();
subplot2(2,1,1,1);
imagesc(vbc{grps{g}},sq(mean(mmJpdf(2,3,:,:,1,g,:),7))');axis('xy');
subplot2(2,1,2,1);
imagesc(vbc{grps{g}},sq(mean(mmJpdf(2,3,:,:,2,g,:),7))');axis('xy');
subplot2(2,1,1,2);
imagesc(vbc{grps{g}},sq(mean(smplCnt(:,:,1,g,:),7))');axis('xy');


figure();
subplot2(2,1,1,1);
plot(sq((mmJpdf(2,3,:,5,1,g,:))));
subplot2(2,1,2,1);
plot(sq((mmJpdf(2,3,:,5,2,g,:))));


figure,
subplot(211);    hold('on');
e = 1;
scl = 10;
yc = numel(vbc{grps{g}(2)});
cmap = cool(yc);
for y = 1:yc
    errorbar(vbc{grps{g}(1)},sq(mean(mmJpdf(2,3,:,y,e,g,:),7)),sq(std(mmJpdf(2,3,:,y,e,g,:),[],7)./scl),'Color',cmap(y,:));
end
subplot(212);    hold('on');
e = 2;
for y = 1:yc
    errorbar(vbc{grps{g}(1)},sq(mean(mmJpdf(2,3,:,y,e,g,:),7)),sq(std(mmJpdf(2,3,:,y,e,g,:),[],7)./scl),'Color',cmap(y,:));
end

%%%>>>



%%%<<< mirrored bin analysis : ABS angvel and hbang

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
sts = 4;
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
ForAllSubplots('xlim([-250,250]);');
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

%%%>>>



%% ---------------------------------------------------------------------------------------


hvangBinEdges = [-0.1,-0.05, -0.03, -0.018 ,-0.009,  0.009, 0.018, 0.03, 0.05, 0.1];
hvangBinEdges = [-0.15,-0.035,  -0.015 , 0.015, 0.035, 0.15];
hvangBinCenters = mean([hvangBinEdges(2:end); hvangBinEdges(1:end-1)]);
edy = numel(hvangBinCenters);
hvangBinInd = discretize(-dhvang,hvangBinEdges);

hpchBinEdges = linspace(-1.4,0,8);
hpchBinCenters = mean([hpchBinEdges(2:end); hpchBinEdges(1:end-1)]);
edy = numel(hpchBinCenters);
hpchBinInd = discretize(dfet,hpchBinEdges);

cdhbang = dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]));

eds = 12;
hbangBinEdges = linspace(-1.2,1.2,eds);
hbangBinCenters = mean([hbangBinEdges(2:end); hbangBinEdges(1:end-1)]);
edx = numel(hbangBinCenters);
hbangBinInd = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),hbangBinEdges);
% $$$ 
% $$$ figure,hist([(-(dhbang(ind)+0.3-0.55*double(~ismember(dtind(ind),[3,4,5]))))],linspace(-pi,pi,500));
% $$$ figure,hist2([(-(dhbang(ind)+0.1-0.2*double(~ismember(dtind(ind),[3,4,5])))),dxyvel(ind)], ...
% $$$              linspace(-pi,pi,100),linspace(-0,3,100));



sts = 1;
ind =   logical(dstcm(:,1))                             ...
        & any(dstcm(:,[3,4,5,6]),2)                  ...
        & dpostI                                        ...
        & duincI                                        ...
        & drberrs<100                                   ...        
        & dhdist<350;


%ind = ind & ~ismember(dtind,[3,4,5]);
%        & any(dstcm(:,stsGrps{sts}),2)                  ...


%%%<<< diagnostics decode lat lon
figure,
for p = 1:8,
subplot(8,1,p);
mind = ind&smMask&dphz==phzBinCenters(p);
%mind = ind&smMask&dphz==phzBinCenters(p)&dhbang<-0.3;
%mind = ind&smMask&dphz==phzBinCenters(p)&dhbang>0.3;
hist2([dErrlon(mind),dErrlat(mind)],...
      linspace(-500,500,50),linspace(-500,500,50));
xlim([-300,300]);
ylim([-300,300]);
Lines([],0,'m');
Lines(0,[],'m');
end
%%%>>>
display = true;
for e = 1:2;
    if display,
        figure();
        sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
    end
    
    tind = ind & dxyvel>5;
    tfet = ferr{e}(tind);       tphz = dphz(tind);
    %txi = hbangBinInd(tind);    tyi = hvangBinInd(tind);
    txi = hbangBinInd(tind);    tyi = hpchBinInd(tind);
    out = {};
    for f = 1:edy,    
        for b = 1:edx
            if display,
                axes(sax(f,b));       
                hold('on');        
            end
            
            indb = tyi==f & txi==b; 
            out{f,b} = hist2([tfet(indb),tphz(indb)],ferrorBinEdges{e},phzBins);
            if display,
                imagesc(ferrorBinEdges{e},...
                        [phzBins,phzBins+2*pi],...
                        imgaussfilt(repmat(out{f,b},[1,2]),[2,0.1])');
                imagesc(ferrorBinEdges{e},phzBins+2*pi,imgaussfilt(out{f,b},[2,0.1])');
                axis('tight');        axis('xy');        xlim([-250,250]);        colormap(gca(),'jet');
                Lines([],pi,'k');        Lines(0,[],'k');
                text(sax(f,b),-220,-pi,['hba:',num2str(round(hbangBinCenters(b),2)),' pch:', ...
                                    num2str(round(hpchBinCenters(f),2))],'Color','w','Rotation',90);
                %text(sax(f,b),-220,-pi,['hba:',num2str(round(hbangBinCenters(b),2)),' ...
                % hva:',num2str(round(hvangBinCenters(f),2))],'Color','w','Rotation',90);            
            end

        end
    end


    out = cf(@(o) bsxfun(@rdivide,o,sum(o,'omitnan')), out);
    %n = cell2mat(cf(@(o) find(cumsum(o(:,3))>0.5,1,'first'), out));
    %m = cf(@(o) permute(sum(bsxfun(@times,o,ferrorBinCenters{1}'),'omitnan'),[1,3,2]), out);
    %m = reshape(cat(1,m{:}),[7,9,8]);


    n = nan([edy,edx,8]);
    if display,
        figure();
    end

    for p = 1:8,
        if display,    subplot(8,1,p);    end;
        for f = 1:edy,
            for b = 1:edx,
                [x,index] = unique(cumsum(out{f,b}(:,phzOrder(p)),'omitnan'));
                index(x==0|x==1) = [];
                x(x==0|x==1) = [];
                [l,li] = unique(cumsum(out{f,b}(:,phzOrder(p)),'omitnan'),'last');
                li(l==0|l==1) = [];
                l(l==0|l==1) = [];
                
                try,
                    n(f,b,p) = (interp1(x,ferrorBinCenters{e}(index),0.5) ...
                                + interp1(l,ferrorBinCenters{e}(li),0.5))./2;
                end
            end
        end
        
        if display,
            imagesc(n(:,:,p));
            colormap('jet');
            if e==1
                caxis([-20,100]);
            else
                caxis([-70,70]);
            end
        end
    end

    if e ==1
        lon=n;
    else
        lat=n;
    end


end

rlat = diff(lat(:,:,[6,2]),1,3);
rlon = diff(lon(:,:,[6,2]),1,3);

rlat = lat(:,:,2);
rlon = lon(:,:,2);


figure();
subplot(1,3,1);
imagesc(sqrt(rlat.^2+rlon.^2))
caxis([0,200]);
subplot(1,3,2);
imagesc(rlat)
caxis([-100,100]);
subplot(1,3,3);
imagesc(rlon);
caxis([0,200]);

figure();
plot(hpchBinCenters,median(sqrt(rlat.^2+rlon.^2)),'-+');

figure
plot(cos(hpchBinCenters).*150,median(sqrt(rlat.^2+rlon.^2),2));
hold('on');
line([50,120],[50,120])


%%%<<< Ego centric phase precession bootstrap (hbang,hvang) 

numIter = 100;
pxavn = zeros( [numel( errorBinCenters ), ...
                numel( phzBinCenters   ), ...
                edx,                      ...
                edy,                      ...
                numel(stsGrps),           ...
                2,                        ...
                numIter,                  ...               
             ]                           ...
           );
qtls = [0.25,0.50,0.75];
qntlavr = nan(  [numel(qtls), ...
                numel( phzBinCenters   ), ...
                edx,                      ...
                edy,                      ...
                numel(stsGrps),           ...
                2,                        ...
                numIter,                  ...               
             ]                           ...
           );
qntlavl = nan(  [numel(qtls), ...
                numel( phzBinCenters   ), ...
                edx,                      ...
                edy,                      ...
                numel(stsGrps),           ...
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
            & duincI                                        ...
            & dhdist<380;
    disp(['[Info] State: ',stsLbls{sts}]);
    
    tdphz = dphz(ind);
    thvangBinInd = hvangBinInd(ind);
    thbangBinInd = hbangBinInd(ind);
    tsmMask = smMask(ind);

    for e = 1:2
        tferr = ferr{e}(ind);
        disp(['[Info] e: ',num2str(e)]);
        for f = 1:edy,
            for b = 1:edx,
                tic
                for i = 1:numIter,
                    indb = thvangBinInd == f & thbangBinInd == b & circshift(tsmMask, randsample(shifts,1));
                    indb(indb) = randn([sum(indb),1]) > 0;
                    pxavn(:,:,b,f,sts,e,i) =               ...
                        histcounts2(tferr(indb),          ...
                                    tdphz(indb),             ...
                                    ferrorBinEdges{e},       ...
                                    phzBins);
                    for p = 1:numel(phzBinCenters),
                        try
                            [x, index] = unique(cumsum(pxavn(:,p,b,f,sts,e,i)./sum(pxavn(:,p,b,f,sts,e,i))),'first');
                            index(x==0|x==1) = [];
                            x(x==0|x==1) = [];
                            qntlavr(:,p,b,f,sts,e,i) = interp1(x,ferrorBinCenters{e}(index),qtls);
                        end
                        try
                            [x, index] = unique(cumsum(pxavn(:,p,b,f,sts,e,i)./sum(pxavn(:,p,b,f,sts,e,i))),'last');
                            index(x==0|x==1) = [];
                            x(x==0|x==1) = [];
                            qntlavl(:,p,b,f,sts,e,i) = interp1(x,ferrorBinCenters{e}(index),qtls);
                        end
                    end
                end
                toc
            end
        end
    end
end
qntlav = cat(8,qntlavr,qntlavl);
qntlav = mean(qntlav,8);

%%%>>>

%%%<<< Ego centric phase precession bootstrap randomized head body angle  (hbang,hvang) 

numIter = 100;
pxavnShuff = zeros( [numel( errorBinCenters ), ...
                numel( phzBinCenters   ), ...
                edx,                      ...
                edy,                      ...
                numel(stsGrps),           ...
                2,                        ...
                numIter,                  ...               
             ]                           ...
           );
qtls = [0.25,0.50,0.75];
qntlavrShuff = nan(  [numel(qtls), ...
                numel( phzBinCenters   ), ...
                edx,                      ...
                edy,                      ...
                numel(stsGrps),           ...
                2,                        ...
                numIter,                  ...               
             ]                           ...
           );
qntlavlShuff = nan(  [numel(qtls), ...
                numel( phzBinCenters   ), ...
                edx,                      ...
                edy,                      ...
                numel(stsGrps),           ...
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
            & duincI                                        ...
            & dhdist<380;
    disp(['[Info] State: ',stsLbls{sts}]);
    
    tdphz = dphz(ind);
    thvangBinInd = hvangBinInd(ind);
    thbangBinInd = hbangBinInd(ind);
    tsmMask = smMask(ind);

    for e = 1:2
        tferr = ferr{e}(ind);
        disp(['[Info] e: ',num2str(e)]);
        for f = 1:edy,
            for b = 1:edx,
                tic
                for i = 1:numIter,
                    cmask = circshift(tsmMask, randsample(shifts,1));
% GET index of mirrored head-body angles.
                    indb =   thvangBinInd == f                                                   ...
                           & ismember(thbangBinInd,[b,edx+1-b])                                  ...
                           & cmask;
% RESAMPLE from population from mirrored head-body angles.
                    indm = false([sum(double(indb)),1]);
                    indm(randsample(sum(double(indb)),...
                                    sum(double(thvangBinInd == f ...
                                               & thbangBinInd == b ...
                                               & cmask)))) = true;
                    indb(indb) = indm;
                    indb(indb) = randn([sum(indb),1]) > 0;
% COMPUTE the JPDF of 
                    pxavnShuff(:,:,b,f,sts,e,i) =                                                     ...
                        histcounts2(tferr(indb),                                                 ...
                                    tdphz(indb),                                                 ...
                                    ferrorBinEdges{e},                                           ...
                                    phzBins);
                    for p = 1:numel(phzBinCenters),
                        try
                            [x, index] = unique(cumsum(pxavnShuff(:,p,b,f,sts,e,i)./sum(pxavnShuff(:,p,b,f,sts,e,i))),'first');
                            index(x==0|x==1) = [];
                            x(x==0|x==1) = [];
                            qntlavrShuff(:,p,b,f,sts,e,i) = interp1(x,ferrorBinCenters{e}(index),qtls);
                        end
                        try
                            [x, index] = unique(cumsum(pxavnShuff(:,p,b,f,sts,e,i)./sum(pxavnShuff(:,p,b,f,sts,e,i))),'last');
                            index(x==0|x==1) = [];
                            x(x==0|x==1) = [];
                            qntlavlShuff(:,p,b,f,sts,e,i) = interp1(x,ferrorBinCenters{e}(index),qtls);
                        end
                    end
                end
                toc
            end
        end
    end
end
qntlavShuff = cat(8,qntlavlShuff,qntlavrShuff);
qntlavShuff = mean(qntlavShuff,8);

%%%>>>


%%%<<< DATA CHECK POINT : qntlav

save(fullfile(MTA_PROJECT_PATH,'analysis','MjgER2016_req20191104_qntlav.mat'),...
     'qntlav',                  ...
     'qntlavShuff',             ...
     'hbangBinCenters',         ...
     'hvangBinCenters',         ...
     'edx','edy',               ...
     'stsGrps',                 ...
     'stsLbls'                  ...
);
load(fullfile(MTA_PROJECT_PATH,'analysis','MjgER2016_req20191104_qntlav.mat'));

%%%>>>


%%%<<< Plot median decoded position (egocentric) bootstraps for thetaPhase([90,270])

ledLoc = {'southwest','northeast'}; 
hba = 1:edx
sts = 1
figure();
for p = [3];    for hva = 1:edy;
    subplot(1,edy,hva);
    hold('on');
    %plot( sq(mean(qntlav( 2, p,
    %hba,hva,sts,1,:),7,'omitnan')),sq(mean(qntlav(2,p,hba,hva,sts,2,:),7,'omitnan')))
    plot( sq(qntlav( 2, 7, hba,hva,sts,1,:)),sq(qntlav(2,7,hba,hva,sts,2,:)),'b')
    plot( sq(qntlav( 2, p, hba,hva,sts,1,:)),sq(qntlav(2,p,hba,hva,sts,2,:)),'c')
    xlim([-120,120]);
    ylim([-120,120]);
end;end;

%%%>>>

%%%<<< Plot bootstrapped median of the median decoded (egocentric) position over thetaPhase([80,130,180,230,270])

cmap = jet(12);
ledLoc = {'southwest','northeast'}; 
hba = 1:edx;
sts = 1;
figure();
for p = [2:7];    for hva = 1:edy;
    subplot(1,edy,hva);
    hold('on');
    plot( sq(mean(qntlav( 2, phzOrder(p),hba,hva,sts,1,:),7,'omitnan')),...
          sq(mean(qntlav( 2, phzOrder(p),hba,hva,sts,2,:),7,'omitnan')),'Color',cmap(p,:));
    xlim([-20,120]);
    ylim([-80,80]);
end;end;

%%%>>>


cmap = jet(5);
ledLoc = {'southwest','northeast'}; 
hba = 2:edx-1;
sts = 1;
figure();
for p = [6];    
    for sts = 1:2,
    for hva = 1:edy;
    subplot(1,edy,hva);
    hold('on');
    plot( sq(qntlav( 2, phzOrder(p),hba,hva,sts,1,:))+randn([numel(hba),numIter])./10,...
          sq(qntlav( 2, phzOrder(p),hba,hva,sts,2,:))+randn([numel(hba),numIter])./10,...
          '.','Color',cmap(sts,:));
    plot( sq(mean(qntlav( 2, phzOrder(p),hba,hva,sts,1,:),7,'omitnan')),...
          sq(mean(qntlav( 2, phzOrder(p),hba,hva,sts,2,:),7,'omitnan')),...
          '-+','Color',cmap(sts,:),'LineWidth',5);
    xlim([-20,120]);
    ylim([-80,80]);
    end
end;end;
legend(stsLbls)

figure()
    plot( sq(mean(qntlav( 2, phzOrder(p),hba,hva,sts,1,:),7,'omitnan')),...
          sq(mean(qntlav( 2, phzOrder(p),hba,hva,sts,2,:),7,'omitnan')),'Color',cmap(sts,:));



sts = 1;
[tt,rr] = cart2pol(diff(qntlav( 2, [7,3],hba,:,sts,1,:),1,2),diff(qntlav( 2, [7,3],hba,:,sts,2,:),1,2));
figure,
subplot(121);
imagesc(circ_mean(sq(tt),[],3)');
subplot(122);
imagesc(circ_std(sq(tt),[],[],3)');

figure,
subplot(121);
imagesc(mean(sq(rr),3)');
subplot(122);
imagesc(std(sq(rr),[],3)');




figure();
plot(hvangBinCenters,circ_mean(circ_mean(sq(tt),[],3)),'.');



figure,
plot(repmat(hbangBinCenters(2:end-1)',[1,100])+randn([numel(hbangBinCenters(2:end-1)),100])./20,sq(mean(sq(tt),2)),'.b');
figure();
boxplot(sq(mean(sq(tt),2))','labels',round(hbangBinCenters(2:end-1),2),'plotstyle','compact');
figure
errorbar(round(hbangBinCenters(2:end-1),2),...
         mean(sq(mean(sq(tt),2)),2),...
         std(sq(mean(sq(tt),2)),[],2).*1.96);




%%%<<< SUPFIG head-body yaw (dhbang) distirbutions across states
[hfig,fig,fax,sax] = set_figure_layout(figure(666061),'A4','portrait',[],1.5,1.5,0.1,0.1);
    
    nx = 3
    ny = 5
    x = 1;
    subplot2(ny,nx,1,x);
        out = histc(-cdhbang(ind&smMask),linspace(-pi,pi,100));
        bar(linspace(-pi,pi,100),out./sum(out),'histc')
        xlim([-2,2]);
        ylim([0,0.1]);
        Lines(hbangBinEdges([2:end-1]),[],'r');        
        title({'Head-Body angle (yaw)','all'});
    for s = 3:6,    
        subplot2(ny,nx,s-1,x);        
            out = histc(-cdhbang(ind&smMask&dstcm(:,s)==s),linspace(-pi,pi,100));
            bar(linspace(-pi,pi,100),out./sum(out),'histc')
            xlim([-2,2]);
            ylim([0,0.1]);
            Lines(hbangBinEdges([2:end-1]),[],'r');        
            title({states{s}});
    end
    xlabel({'radians','(0.06rad/bin)'});
    ylabel('probability');

% SUPFIG head-body yaw angvel (dhvang) distirbutions across states    
    x = 2;
    subplot2(ny,nx,1,x);
        out = histc(-dhvang(ind&smMask),linspace(-0.3,0.3,100));
        bar(linspace(-0.3,0.3,100),out./sum(out),'histc')
        xlim([-0.3,0.3]);
        ylim([0,0.2]);
        Lines(hvangBinEdges,[],'r');        
        title({'Head-Body angvel (d(yaw)/dt)','all'});
    for s = 3:6,    
        subplot2(ny,nx,s-1,x);
            out = histc(-dhvang(ind&smMask&dstcm(:,s)==s),linspace(-0.5,0.5,100));
            bar(linspace(-0.3,0.3,100),out./sum(out),'histc')
            xlim([-0.3,0.3]);
            ylim([0,0.2]);
            Lines(hvangBinEdges,[],'r');        
            title({states{s}});
    end
    xlabel({'radians/sample(@120Hz)','(0.06rad/sample/bin)'});
    ylabel('probability');    

% SUPFIG JPDF head-body yaw vs angular velocity distirbutions across states
    x = 3;
    subplot2(ny,nx,1,x);
        out = hist2([-cdhbang(ind&smMask),dhvang(ind&smMask)],...
                    linspace(-pi,pi,50),linspace(-0.3,0.3,50));
        imagesc(linspace(-pi,pi,50),linspace(-0.3,0.3,50),out'./sum(out(:)));
        Lines(hbangBinEdges(2:end-1),[],'r');
        Lines([],hvangBinEdges(2:end-1),'r');
        ylim([-0.15,0.15]);        
        xlim([-1.5,1.5]);
        axis('xy');
        caxis([0,0.015]);
        title({'Head-Body yaw VS','Head-Body angvel (d(yaw)/dt)','all'});
    for s = 3:6,    
        subplot2(ny,nx,s-1,x);
            out = hist2([-cdhbang(ind&smMask&dstcm(:,s)==s),dhvang(ind&smMask&dstcm(:,s)==s)],...
                        linspace(-pi,pi,50),linspace(-0.3,0.3,50));
            imagesc(linspace(-pi,pi,50),linspace(-0.3,0.3,50),out'./sum(out(:)));
            ylim([-0.15,0.15]);        
            xlim([-1.5,1.5]);
            Lines(hbangBinEdges(2:end-1),[],'r');
            Lines([],hvangBinEdges(2:end-1),'r');
            axis('xy');
            caxis([0,0.015]);
            title({states{s}});
    end
    xlabel({'radians','(0.06rad/bin)'});
    ylabel('rad/sample');
    sax =  gca();    
    cax = colorbar(sax);
    drawnow();
    cax.Position(1) = sum(sax.Position([1,3]));
    drawnow();    
    ylabel(cax,'Probability');
    % Because matlab sometimes sucks
    cax.Position(1) = sum(sax.Position([1,3]));        

figdir = create_directory(fullfile(MTA_PROJECT_REPORT_PATH,'decode/phase_precession'));
%print(hfig, '-depsc2', fullfile(figdir,'subfig_headBodyYaw_distributions.eps'));
print(hfig, '-dpng', fullfile(figdir,'subfig_headBodyYaw_distributions.png'));
    
%%%>>>
    
    

    
    
figure();
p = 2;
sts = 1;
for hva = 1:5,
    subplot(1,5,hva);
    boxplot(sq(qntlav(2,phzOrder(p),hba,hva,sts,2,:))','labels',round(hbangBinCenters(2:end-1),2),'plotstyle','compact');
end

figure
errorbar(round(hbangBinCenters(2:end-1),2),...
         sq(mean(qntlav(2,phzOrder(p),hba,hva,sts,2,:),7,'omitnan')),...
         sq(std(qntlav(2,phzOrder(p),hba,hva,sts,2,:),[],7,'omitnan')).*1.96);
            

figure,
plot(linspace(-pi/2,mean(circ_mean(sq(tt),[],3),2),'.');


[R,P] = corrcoef(hbangBinCenters(2:end-1),mean(circ_mean(sq(tt),[],3),2));
[R,P] = corrcoef(hbangBinCenters(2:end-1),mean(circ_mean(sq(tt),[],3),2));


mtt = circ_mean(sq(tt),[],3);
figure();
plot(mean(mtt)','.');


cmap = cool(5);
figure();
hold('on');
for j = 1:5,
    plot(hbangBinCenters(2:end-1),mtt(:,j),'.','Color',cmap(j,:),'MarkerSize',10);
end
line([-0.8,1],[-0.8,1])



cmap = cool(5);
figure();hold('on');
for j = 1:5;
    plot(diff(circ_mean(sq(tt(:,:,:,j,:,:,:)),[],2)),'Color',cmap(j,:));
end
plot(mean(diff(circ_mean(sq(tt),[],3))'),'k--','LineWidth',2)



figure
p = 6;
hba = 2:edx-1;
hva = 1:edy;
for sts = 1:5
    subplot2(2,5,1,sts);
    imagesc(round(circ_rad2ang(hbangBinCenters(hba))),...
            hvangBinCenters(hva),...
            sq(mean(qntlav( 2, p, hba,hva,sts,1,:),7,'omitnan'))'./10);
    xlabel('Head-Body angle');
    ylabel('Head-Body angvel');    
    cax = colorbar();
    ylabel(cax,'cm');    
    colormap('cool');
    caxis([-1,12]);
    title(stsLbls{sts});
    subplot2(2,5,2,sts);
    imagesc(round(circ_rad2ang(hbangBinCenters(hba))),...
            hvangBinCenters(hva),...
            sq(mean(qntlav( 2, p, hba,hva,sts,2,:),7,'omitnan'))'./10);
    xlabel('Head-Body angle');
    ylabel('Head-Body angvel');    
    cax = colorbar();
    ylabel(cax,'cm');
    colormap('cool');
    caxis([-7,7]);
end

figure();
sax = reshape(tight_subplot(8,5,[0.01,0.01],[0.15,0.15],[0.15,0.2]),[5,8])';
fax = axes('Position',[0,0,1,1],'Visible','off');
hba = 2:edx-1;
hva = 1:edy;
for p = 1:8;
    for sts = 1:5
        axes(sax(p,sts));           
        qNDiff = reshape(bsxfun(@minus,...
                                permute(qntlav( 2, phzOrder(p), hba,hva,sts,2,:),[1,2,3,4,5,7,6]), ...
                                qntlav( 2, phzOrder(p), fliplr(hba),hva,sts,2,:)),...
                         numel(hba),numel(hva),[]);
        qSDiff = reshape(bsxfun(@minus,...
                                permute(qntlavShuff( 2, phzOrder(p), hba,hva,sts,2,:),[1,2,3,4,5,7,6]), ...
                                qntlavShuff( 2, phzOrder(p), fliplr(hba),hva,sts,2,:)),...
                         numel(hba),numel(hva),[]);
        
        zscr = (mean(qNDiff,3)-mean(qSDiff,3))./std(qSDiff,[],3);;
        %pval = 2*normcdf(zscr);

        imagesc(round(circ_rad2ang(hbangBinCenters(hba))),...
                hvangBinCenters(hva),...
                sq(mean(qntlav( 2, phzOrder(p), hba,hva,sts,2,:),7,'omitnan'))'./10);        
        %zscr');

        caxis([-8,8]);
        colormap('jet');    
        if p == 1,
            title(stsLbls{sts});
        end
        if sts ==1 && p == 8,
            ylabel({'Head-Body','angvel'});    
        end
        if sts ==3 && p == 8,        
            xlabel({'Head-Body','angle'});
        end
        
        
        if p ~= 8,
            set(gca(),'XTickLabel',{});
        end
        if sts ~= 1
            set(gca(),'YTickLabel',{});
        end
        if sts==5         
            if p==8,
                cax = colorbar();
                %ylabel(cax,'cm');
                ylabel(cax,'Z-Score');
                cax.Position(1) = [sum(sax(p-1,sts).Position([1,3])+0.02)];
            end
            pause(0.5);
            axes(fax);
            text(sum(sax(p,sts).Position([1,3]))+0.01,...
                 sum(sax(p,sts).Position([2,4]).*[1,0.5]),...
                 num2str(round(circ_rad2ang(phzBinCenters(phzOrder(p))+2.*pi.*double(phzBinCenters(phzOrder(p))<0)))),...
                 'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90);
        end
    end
end
axes(fax);
text(0.85,0.5,'Theta Phase','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90);
text(0.5,0.925,'Lateral Phase Precession','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',0);




figure();
sax = reshape(tight_subplot(8,5,[0.01,0.01],[0.15,0.15],[0.15,0.2]),[5,8])';
fax = axes('Position',[0,0,1,1],'Visible','off');
hba = 2:edx-1;
hva = 1:edy;
for p = 1:8;
    for sts = 1:5
        axes(sax(p,sts));           
        qNDiff = reshape(bsxfun(@minus,...
                                permute(qntlav( 2, phzOrder(p), hba,hva,sts,1,:),[1,2,3,4,5,7,6]), ...
                                qntlav( 2, phzOrder(p), fliplr(hba),hva,sts,1,:)),...
                         numel(hba),numel(hva),[]);
        qSDiff = reshape(bsxfun(@minus,...
                                permute(qntlavShuff( 2, phzOrder(p), hba,hva,sts,1,:),[1,2,3,4,5,7,6]), ...
                                qntlavShuff( 2, phzOrder(p), fliplr(hba),hva,sts,1,:)),...
                         numel(hba),numel(hva),[]);
        
        zscr = (mean(qNDiff,3)-mean(qSDiff,3))./std(qSDiff,[],3);;
        %pval = 2*normcdf(zscr);

% $$$         imagesc(round(circ_rad2ang(hbangBinCenters(hba))),...
% $$$                 hvangBinCenters(hva),...
% $$$                 zscr');
        
        imagesc(round(circ_rad2ang(hbangBinCenters(hba))),...
                hvangBinCenters(hva),...
                sq(mean(qntlav( 2, phzOrder(p), hba,hva,sts,1,:),7,'omitnan'))'./10);
        %sq(mean(qntlav( 2, phzOrder(p), hba,hva,sts,1,:),7,'omitnan'))'./10);        
        caxis([-5,10]);
        colormap('jet');    
        if p == 1,
            title(stsLbls{sts});
        end
        if sts ==1 && p == 8,
            ylabel({'Head-Body','angvel'});    
        end
        if sts ==3 && p == 8,        
            xlabel({'Head-Body','angle'});
        end
        
        
        if p ~= 8,
            set(gca(),'XTickLabel',{});
        end
        if sts ~= 1
            set(gca(),'YTickLabel',{});
        end
        if sts==5         
            if p==8,
                cax = colorbar();
                %ylabel(cax,'cm');
                ylabel(cax,'Z-Score');
                cax.Position(1) = [sum(sax(p-1,sts).Position([1,3])+0.02)];
            end
            pause(0.5);
            axes(fax);
            text(sum(sax(p,sts).Position([1,3]))+0.01,...
                 sum(sax(p,sts).Position([2,4]).*[1,0.5]),...
                 num2str(round(circ_rad2ang(phzBinCenters(phzOrder(p))+2.*pi.*double(phzBinCenters(phzOrder(p))<0)))),...
                 'VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90);
        end
    end
end
axes(fax);
text(0.85,0.5,'Theta Phase','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',90);
text(0.5,0.925,'Forward Phase Precession','VerticalAlignment','middle','HorizontalAlignment','center','Rotation',0);







figure
for sts = 1:5
    subplot2(2,5,1,sts);
    imagesc(sq(mean(diff(qntlav( 2, [7,4], :,:,sts,1,:),1,2),7,'omitnan'))');
    colorbar();colormap('cool');caxis([-30,120]);title(stsLbls{sts});
    subplot2(2,5,2,sts);
    imagesc(sq(mean(diff(qntlav( 2, [7,4], :,:,sts,2,:),1,2),7,'omitnan'))');
    colorbar();colormap('cool');caxis([-70,70]);
end


figure,imagesc(sq(mean(qntlav( 2, 7, :,1:5,1,2,:),7,'omitnan'))');colormap cool
figure,imagesc(sq(mean(qntlav( 2, 7, :,1:5,1,1,:),7,'omitnan'))');colormap cool
plot( sq(mean(qntlav( 2, p, hba,hva,1,1,:),7,'omitnan')),sq(mean(qntlav(2,p,hba,hva,1,2,:),7,'omitnan')))


%%%<<< Compute z-score 

% does the sign of the hba matter?
% bootstrapped median decoded lateral displacement relative to the head
%    compare the sign of the head-body angle for each group of displacements
%    compare the sign of the head-body angular velocity for each group of velocities
hba = 1:numel(hbangBinCenters);

randsample(qntlav(:,p,hba,hav,sts,fet,itr),50;

%%%>>>


%% STARTFIG ------------------------------------------------------------------------------

[hfig,fig,fax,sax] = set_figure_layout(figure(666008),'A4','portrait',[],1.5,1.5,0.1,0.1);


ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[9:11])),2)            ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<320;


vLbl = {};
vEds = {};
vCtr = {};
vInd = {};
vBinLbl = {};
varBinCnt = 3;

vLbl{end+1} = 'hba';
vEds{end+1} = [-1.2,-0.3,0.3,1.2];  % length: varBinCnt
%vEds{end+1} = linspace(0,1.2,6);
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),vEds{end});
%vBinLbl{end+1} = {'L','CL','C','CR','R'};
vBinLbl{end+1} = {'Left','Center','Right'};

vLbl{end+1} = 'hvl';
vEds{end+1} = [-80,-5,5,80]; % length: varBinCnt
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(dfhrvl,vEds{end});
vBinLbl{end+1} = {'CCW','0','CW'};

vLbl{end+1} = 'hvf';
vEds{end+1} = [-10,10,25,80]; % length: varBinCnt
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(dfhrvf,vEds{end});
vBinLbl{end+1} = {'rest','slow','fast'};
% $$$ 
% $$$ vLbl{end+1} = 'hbp';
% $$$ vEds{end+1} = [-1.5,-0.75,-0.25,0.5]; % length: varBinCnt
% $$$ vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
% $$$ vInd{end+1} = discretize(dpch,vEds{end});
% $$$ vBinLbl{end+1} = {'down','level','up'};


% COMPUTE JPDF(FERROR, PHZ | HBA, HRVL)
tdphz = dphz(ind);

out =  zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters),varBinCnt, varBinCnt, 2,numel(vLbl),numel(vLbl)]);
for v = 1;%:numel(vLbl),
    for b = 1:numel(vLbl),
        if v==b, continue, end        
% SETUP conditioned vars
        edx = numel(vCtr{v});
        edy = numel(vCtr{b});
        txi = vInd{v}(ind);
        tyi = vInd{b}(ind);

        for e = 1:2;
% GET dependent var { ego-centric phase precession: forward, lateral }
            tferr = ferr{e}(ind);
            if e==2, tferr = tferr+8; end
% COMPUTE conditional means
            for x = 1:edx
                for y = 1:edy,            
                    indb = tyi==y & txi==x;
                    out(:,:,x,y,e,v,b) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
                end
            end
        end        
    end
end


mout =  zeros([1,numel(phzBinCenters),varBinCnt, varBinCnt, 2,numel(vLbl),numel(vLbl)]);
for v = 1,%:numel(vLbl),
    for b = 1:numel(vLbl),
        if v==b, continue, end
        for y = 1:edy,    
            for x = 1:edx
                for p = 1:numel(phzBinCenters),
                   for e = 1:2,
                        try
                            [g, index] = unique(cumsum(out(:,p,x,y,e,v,b)./sum(out(:,p,x,y,e,v,b))),'last'); 
                            index(g==0|g==1) = [];
                            g(g==0|g==1) = [];
                            mout(:,p,x,y,e,v,b) = interp1(g,ferrorBinCenters{e}(index),0.5);
                        end
                    end
                end
            end
        end
    end 
end

% out[ prj, phz, x, y, epp, var1, var2 ]


[hfig,fig,fax,sax] = set_figure_layout(figure(666010),'A4','portrait',[],1.25,1.25,0.1,0.1);

%       hba,hvl  hvl,hvf  hvf,hbp  hba,hvf
vGrp = { [1,2],   [3,2],   [3,1]};
vGrp = { [1,2],   [1,3]};
cmap = 'bcr';
cmap = 'bgr';    
for g = 1:numel(vGrp);
    v = vGrp{g}(1); 
    b = vGrp{g}(2); 
    % PLOT JPDF(FERROR, PHZ | HBA, HRVL)
    [yind, yOffset, xind, xOffset] = deal(1,0, 1, 0);
    for e = 1:2,
        for x = 1:edx
            for y = 1:edy,                
                [yind, yOffSet, xind, xOffSet] = deal(y,                                        ...
                                                      y*fig.subplot.height/2-(g-1)*4,           ...
                                                      x+3*double(e==2),                         ...
                                                      fig.subplot.width/3+double(e==2)./4-1);
                sax(end+1) = axes('Units','centimeters',                                        ...
                                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                    fig.page.ypos(yind)+yOffSet,                      ...
                                    fig.subplot.width ,                               ...
                                    fig.subplot.height/2],                            ...
                                  'FontSize', 8,                                                ...
                                  'LineWidth',1);
                hold(sax(end),'on');
                imagesc(ferrorBinEdges{e},                                           ...
                        [phzBins,phzBins+2*pi],                                      ...
                        imgaussfilt(repmat(out(:,:,x,4-y,e,v,b),[1,2]),[3,0.2])');
                axis('tight');  axis('xy');  ylim([0,2.*pi]);
                if e == 1,  xlim([-200,300]);  else  xlim([-250,250]);  end
                colormap(sax(end),'jet');
                Lines([],pi,'k');
                Lines(0,[],'k');
                sax(end).XTickLabel = {};
                sax(end).YTickLabel = {};
                if y == 1 && x == 2 && e == 1 && g==1,
                    title('Forward');
                elseif y == 1 && x == 2 && e == 2 && g==1,
                    title('Lateral');
                end
                if y==3,
                    xlabel(sax(end),vBinLbl{v}(x));
                end
                if x==1&&y==2&&e==1,
                    ylabel(sax(end),[vLbl{vGrp{g}(1)},' vs ',vLbl{vGrp{g}(2)}]);
                end
                if x==3 & e==2,
                    sax(end).YAxisLocation = 'right';
                    yh = ylabel(sax(end),vBinLbl{b}(4-y),'Rotation',0,'Color',cmap(4-y));
                    yh.Position(1:2) = [yh.Position(1)+220,6];
                    %uistack(yh,'top');
                end
            end
        end
    end


    [yind, yOffSet, xind, xOffSet] = deal(y, ...
                                          y*fig.subplot.height/2-(g-1)*4,...
                                          x+4*double(e==2), ...
                                          fig.subplot.width/2+double(e==2));
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                        fig.page.ypos(yind)+yOffSet,                      ...
                        fig.subplot.width*3+fig.subplot.horizontalPadding.*2,   ...
                        fig.subplot.height/2*3+2*fig.subplot.verticalPadding],...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');

    offset = [-50:50:50];
    for j = 1:varBinCnt,
        for i = 1:varBinCnt,
            plot(sq(mout(1,phzOrder(fliplr([2,6])),j,i,2,v,b))+offset(j),...
                 sq(mout(1,phzOrder(fliplr([2,6])),j,i,1,v,b)),[cmap(i),'-'],'LineWidth',1.5)
            plot(sq(mout(1,phzOrder(fliplr(6)),  j,i,2,v,b))+offset(j),...
                 sq(mout(1,phzOrder(fliplr(6))  ,j,i,1,v,b)),[cmap(i),'o'],'LineWidth',2)
        end
    end
    xlim([offset(1)-diff(offset(1:2)),offset(end)+diff(offset(1:2))]);
    ylim([-50,125])
    Lines(-offset,[],'k');
    Lines([],0,'k');
    Lines(offset,[],'k');
    Lines(0,[],'k');
    sax(end).XTick = offset;
    sax(end).XTickLabel = vBinLbl{v};
    sax(end).YAxisLocation = 'right';
end




axes(sax(end));
yOffSet = -0.5;
xOffSet = fig.subplot.width/2;
for s = 1:numel(stsLbls),
    [yind, yOffSet, xind, xOffSet] = deal(4,yOffSet,s,xOffSet);
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width ,                               ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    imagesc(sq(mean(qntlav( 2, 3, 2:end-1,1:5,s,1,:),7,'omitnan'))');
    sax(end).XTickLabel = {};
    sax(end).YTickLabel = {};
    colormap(sax(end),'jet');caxis([-30,120]);title(stsLbls{s});axis('tight');axis('xy');
end

yOffSet = -0.5;
xOffSet = fig.subplot.width/2;
for s = 1:numel(stsLbls),
    [yind, yOffSet, xind, xOffSet] = deal(5,yOffSet,s,xOffSet);
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                                  fig.page.ypos(yind)+yOffSet,                      ...
                                  fig.subplot.width ,                               ...
                                  fig.subplot.height],                              ...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');
    imagesc(sq(mean(qntlav( 2, 3, 2:end-1,1:5,s,2,:),7,'omitnan'))');    
    sax(end).XTickLabel = {};
    sax(end).YTickLabel = {};
    colormap(sax(end),'jet');caxis([-70,70]);title(stsLbls{s});axis('tight');axis('xy');
end












figure();
sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
tind = ind;
tfet = ferr{e}(tind);       tphz = dphz(tind);
tvel = dxyvel(tind);
%txi = hbangBinInd(tind);    tyi = hvangBinInd(tind);
txi = hbangBinInd(tind);    tyi = hpchBinInd(tind);
out = {};
for f = 1:edy,    
    for b = 1:edx
        axes(sax(f,b));       
        hold('on');        
        indb = tyi==f & txi==b; 
        bar(linspace(-1,2,30),histc(log10(tvel(indb)),linspace(-1,2,30)),'histc');
        Lines(0,[],'r');
        Lines(1,[],'r');
        xlim([-1,2]);
    end
end



figure();
sax = reshape(tight_subplot(edy,edx,0.001,0,0),[edx,edy])';
sts = 3
ind =   logical(dstcm(:,1))              ...
        & any(logical(dstcm(:,stsGrps{sts})),2)              ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<360;

tind = ind;
txi = hbangBinInd(tind);    tyi = hvangBinInd(tind);
thvang = dhvang(tind);
out = {};
for f = 1:edy,    
    for b = 1:edx
        axes(sax(f,b));       
        hold('on');        
        indb = tyi==f & txi==b; 
        bar(linspace(0,0.5,100),log10(histc(abs(thvang(indb)),linspace(0,0.5,100))),'histc');
        xlim([0,0.5]);
    end
end

linkaxes(findobj(gcf(),'Type','Axes'),'x');





% GLM fitting of data


sts = 1;
ind =   logical(dstcm(:,1))                             ...
        & any(dstcm(:,[3,4,5,6]),2)                     ...
        & dpostI                                        ...
        & duincI                                        ...
        & drberrs<100                                   ...        
        & dhdist<300                                    ...
        & smMask;

[B,DEV,STATS] = glmfit( [-cdhbang(ind),dhvang(ind),],[dErrlat(ind)]);
glmBinEdges = {[-1.5:0.1:1.5],[-0.15:0.01:0.15]};
glmBinCenters = cf(@(b) mean([b(1:end-1),b(2:end)]), glmBinEdges);
xtest = cell([1,numel(glmBinEdges)]);               
[xtest{:}] = ndgrid(glmBinEdges{:});
xtest = cf(@(x) x(:), xtest);
xtest = cat(2,xtest{:});
yfit = glmval(B,xtest,'identity');
figure,
imagesc(glmBinEdges{:},reshape(yfit,cell2mat(cf(@numel,glmBinEdges)))');
axis xy
caxis([-70,70]);
colormap('jet');
colorbar();


sind = ind&dphz==phzBinCenters(3);
sind = ind&dphz==phzBinCenters(7);
sind(sind) = randn([sum(sind),1])>0;

Y = [dErrlon(sind),dErrlat(sind)];
Xmat = [ones([sum(sind),1]),-dfhrvf(sind),dfhrvl(sind),-cdhbang(sind),dpch(sind)];
d = eye(size(Y,2));
Xcell = cell([1,size(Xmat,1)]);
for i = 1:size(Xmat,1),
    Xcell{i} = [kron([Xmat(i,:)],d)];
end

[beta,sigma,E,V] = mvregress( Xcell,Y);

se = reshape(sqrt(diag(V)),2,5)'

B = reshape(beta,2,5)'

z = E/chol(sigma);
figure()
plot(z(:,1),z(:,2),'.')
title('Standardized Residuals')
hold on

% Overlay standard normal contours
z1 = linspace(-5,5);
z2 = linspace(-5,5);
[zx,zy] = meshgrid(z1,z2);
zgrid = [reshape(zx,100^2,1),reshape(zy,100^2,1)];
zn = reshape(mvnpdf(zgrid),100,100);
[c,h] = contour(zx,zy,zn);
clabel(c,h)



B = {};
DEV = {};
STATS = {};
GLML = {'Lon-Elat'};
[B{end+1},DEV{end+1},STATS{end+1}] = glmfit( [dhrvf(tind)],dErrlat(tind));
GLML = {'LonLat-Elat'};
GLML{end+1} = {'LonYaw-Elat'};
[B{end+1},DEV{end+1},STATS{end+1}] = glmfit( [dhrvf(tind),-cdhbang(tind)],dErrlat(tind));
GLML{end+1} = {'LonLatYaw-Elat'};
[B{end+1},DEV{end+1},STATS{end+1}] = glmfit( [dhrvf(tind),dhrvl(tind),-cdhbang(tind)],dErrlat(tind));
GLML{end+1} = {'Yaw-Elat'};
[B{end+1},DEV{end+1},STATS{end+1}] = glmfit( [-cdhbang(tind)],dErrlat(tind));

STATS = cat(1,STATS{:})


[BR,BINT,R,RINT,STATSR] = regress(dErrlat(tind),[ones([sum(tind),1]),dhrvl(tind)]);
STATSR
[BR,BINT,R,RINT,STATSR] = regress(dErrlat(tind),[ones([sum(tind),1]),dhvang(tind)]);
[BR(2),STATSR]
[BR,BINT,R,RINT,STATSR] = regress(dErrlat(tind),[ones([sum(tind),1]),-cdhbang(tind)]);
[BR(2),STATSR]
[BR,BINT,R,RINT,STATSR] = regress(dErrlat(sind),[ones([sum(sind),1]),-cdhbang(sind)]);
[BR(2),STATSR]


%[B{end+1},DEV{end+1},STATS{end+1}] = glmfit( [dhrvf(sind),dhrvl(sind)],dErrlat(tind)-dErrlat(sind));
% $$$ glmBinEdges = {[-60:3:120],[-70:3:70]};
% $$$ glmBinCenters = cf(@(b) mean([b(1:end-1),b(2:end)]), glmBinEdges);
% $$$ xtest = cell([1,numel(glmBinEdges)]);               
% $$$ [xtest{:}] = ndgrid(glmBinEdges{:});
% $$$ xtest = cf(@(x) x(:), xtest);
% $$$ xtest = cat(2,xtest{:});
% $$$ yfit = glmval(B,xtest,'identity');
% $$$ figure,
% $$$ imagesc(glmBinEdges{:},reshape(yfit,cell2mat(cf(@numel,glmBinEdges)))');
% $$$ axis xy
% $$$ caxis([-70,70]);
% $$$ colormap('jet');
% $$$ colorbar();


B = {};
DEV = {};
STATS = {};
GLML{end+1} = {'LonYaw-Elat'};
[B{end+1},DEV{end+1},STATS{end+1}] = glmfit( [-cdhbang(tind)],dErrlat(tind) );
glmBinEdges = {[-1.5:0.1:1.5]};
glmBinCenters = cf(@(b) mean([b(1:end-1),b(2:end)]), glmBinEdges);
xtest = cell([1,numel(glmBinEdges)]);               
[xtest{:}] = ndgrid(glmBinEdges{:});
xtest = cf(@(x) x(:), xtest);
xtest = cat(2,xtest{:});
yfit = glmval(B{end},xtest,'identity');

figure,hold('on');
%plot(-cdhbang(tind),dErrlat(tind),'.')
hist2([-cdhbang(tind),dErrlat(tind)],linspace(-pi/2,pi/2,30),linspace(-300,300,30));
plot(glmBinEdges{:},yfit,'r');
sts = 1;
p = 2;
sclr = ['bw';'mg';'rc'];
plot(round(hbangBinCenters(2:end-1),2),sq(mean(qntlav(2,phzOrder(p),hba,:,sts,2,:),4)),sclr(sts,1))
plot(round(hbangBinCenters(2:end-1),2),sq(mean(qntlav(1,phzOrder(p),hba,:,sts,2,:),4)),sclr(sts,2))
plot(round(hbangBinCenters(2:end-1),2),sq(mean(qntlav(3,phzOrder(p),hba,:,sts,2,:),4)),sclr(sts,2))
axis('tight');
Lines([],0,'k');Lines([],100,'k');Lines([],-100,'k');






B = {};
DEV = {};
STATS = {};
GLML{end+1} = {'Yaw-Elon'};
[B{end+1},DEV{end+1},STATS{end+1}] = glmfit( [-cdhbang(tind)],dErrlat(tind) );
glmBinEdges = {[-1.5:0.1:1.5]};
glmBinCenters = cf(@(b) mean([b(1:end-1),b(2:end)]), glmBinEdges);
xtest = cell([1,numel(glmBinEdges)]);               
[xtest{:}] = ndgrid(glmBinEdges{:});
xtest = cf(@(x) x(:), xtest);
xtest = cat(2,xtest{:});
yfit = glmval(B{end},xtest,'identity');

figure,hold('on');
%plot(-cdhbang(tind),dErrlat(tind),'.')
hist2([-cdhbang(tind),dErrlat(tind)],linspace(-pi/2,pi/2,30),linspace(-300,300,30));
plot(glmBinEdges{:},yfit,'r');
e = 2;
sts = 2;
p = 2;
sclr = ['bw';'mg';'rc'];
plot(round(hbangBinCenters(2:end-1),2),sq(mean(qntlav(2,phzOrder(p),hba,:,sts,e,:),4)),sclr(sts,1))
plot(round(hbangBinCenters(2:end-1),2),sq(mean(qntlav(1,phzOrder(p),hba,:,sts,e,:),4)),sclr(sts,2))
plot(round(hbangBinCenters(2:end-1),2),sq(mean(qntlav(3,phzOrder(p),hba,:,sts,e,:),4)),sclr(sts,2))
axis('tight');
Lines([],0,'k');Lines([],100,'k');Lines([],-100,'k');




sts = 2;
p = 2;
ind =   logical(dstcm(:,1))                             ...
        & any(dstcm(:,stsGrps{sts}),2)                  ...
        & dpostI                                        ...
        & duincI                                        ...
        & drberrs<100                                   ...        
        & dhdist<350                                    ...
        & smMask;
tind = ind&dphz==phzBinCenters(phzOrder(p));
figure,
for a = 2:10;
    subplot(3,3,a-1);
    hold('on');
    out = hist2([dErrlat(tind&hbangBinInd==a),dErrlon(tind&hbangBinInd==a)],linspace(-250,250,35),linspace(-150,350,35));
    imagesc(linspace(-250,250,35),linspace(-150,350,35),imgaussfilt(out,[2,2])');
    plot(mean(sq(mean(qntlav(2,phzOrder(p),hba,:,sts,2,:),4)),2),...
         mean(sq(mean(qntlav(2,phzOrder(p),hba,:,sts,1,:),4)),2),...
         '-bo','LineWidth',2,'MarkerSize',5,'MarkerFaceColor','g','MarkerEdgeColor','k')
    %plot(sq(mean(qntlav(2,phzOrder(p),a,:,sts,2,:),4)),sq(mean(qntlav(2,phzOrder(p),a,:,sts,1,:),4)),'.g')

    X = sq(mean(qntlav(2,phzOrder(p),a,:,sts,[2,1],:),4));
    Xm = mean(X,2);
    X = bsxfun(@minus,X,Xm);
    error_ellipse(X*X'./(size(X,2)-1),Xm,'conf',0.99,'style','m')
    axis('tight');
    Lines([],100,'k');
    Lines([],0,'k');
    Lines(0,[],'k');    
end







imagesc(glmBinEdges{:},reshape(yfit,cell2mat(cf(@numel,glmBinEdges)))');
axis xy
caxis([-70,70]);
colormap('jet');
colorbar();

GLML{end+1} = {'Yaw-Elat'};
[B{end+1},DEV{end+1},STATS{end+1}] = glmfit( [-cdhbang(ind)],[dErrlat(ind)]);


gind = -cdhbang(sind)<-0.4;
figure,plot(STATS{1}.resid(gind),STATS{3}.resid(gind),'.')

figure();hold('on');
plot(-cdhbang(sind),dErrlat(sind),'.')
figure,hold('on');
plot(-cdhbang(tind),dErrlat(tind),'.')
plot(-cdhbang(tind),STATS{5}.resid,'.')


figure();hold('on');
plot(-cdhbang(tind),STATS{1}.resid-STATS{3}.resid,'.')
plot(linspace(-pi,pi,100),polyval([150/pi,0],linspace(-pi,pi,100)),'-r');

%[B,DEV,STATS] = glmfit( [-cdhbang(ind),dfet(ind)],[dErrlat(ind)]);


%[xtest{:}] = ndgrid([-1:0.1:1],[-1.5:0.01:0],[0:2:60]);
%[xtest{:}] = ndgrid([-1.5:0.01:0],[0:2:60]);
%[xtest{:}] = ndgrid([-1.5:0.1:1.5],[-1.5:0.01:0]);

chi2cdf(DEV{3}-DEV{2},STATS{3}.dfe-STATS{2}.dfe)

figure();
    eds = -800:10:800;
    hist(STATS{1}.resid,eds);
    hold('on');
    hax = bar(eds,histc(STATS{4}.resid,eds),'histc');
    hax.FaceColor = 'r';
    hax.EdgeColor = 'r';
    
    
    
    
%%%<<< MULTIDIMENSIONAL

ind =   logical(dstcm(:,1))                             ...
        & any(logical(dstcm(:,[3,4,5,6])),2)            ...
        & ~any(logical(dstcm(:,[7,8])),2)               ...
        & dpostI                                        ...
        & duincI ...
        & dhdist<380;

vLbl = {};
vEds = {};
vCtr = {};
vInd = {};
vBinLbl = {};
varBinCnt = 3;

vLbl{end+1} = 'hba';
vEds{end+1} = [-1.2,-0.2,0.2,1.2];  % length: varBinCnt
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(-(dhbang+0.2-0.4*double(~ismember(dtind,[3,4,5]))),vEds{end});
vBinLbl{end+1} = {'left','center','right'};

vLbl{end+1} = 'hvl';
vEds{end+1} = [-60,-5,5,60]; % length: varBinCnt
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(dfhrvl,vEds{end});
vBinLbl{end+1} = {'left','0','right'};

vLbl{end+1} = 'hvf';
vEds{end+1} = [-5,5,15,80]; % length: varBinCnt
vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
vInd{end+1} = discretize(dfhrvf,vEds{end});
vBinLbl{end+1} = {'rest','slow','fast'};
% $$$ 
% $$$ vLbl{end+1} = 'hbp';
% $$$ vEds{end+1} = [-1.5,-0.75,-0.25,0.5]; % length: varBinCnt
% $$$ vCtr{end+1} = mean([vEds{end}(2:end); vEds{end}(1:end-1)]);
% $$$ vInd{end+1} = discretize(dpch,vEds{end});
% $$$ vBinLbl{end+1} = {'down','level','up'};


% COMPUTE JPDF(FERROR, PHZ | HBA, HRVL)
tdphz = dphz(ind);

spc = [1,2,3];            
edc = cellfun(@numel,vCtr(spc));
mvjpdf =  zeros([numel(ferrorBinCenters{e}),numel(phzBinCenters), edc, 2]);            

% SETUP conditioned vars
bsi = {};
for s = spc,
    bsi{end+1} = vInd{s}(ind);
end
        

for e = 1:2;
% GET dependent var { ego-centric phase precession: forward, lateral } 
   tferr = ferr{e}(ind);
% COMPUTE conditional means
    for x = 1:edc(1),
        for y = 1:edc(2)
            for z = 1:edc(3),
                indb = bsi{1}==x & bsi{2}==y & bsi{3}==z;
                mvjpdf(:,:,x,y,z,e) = hist2([tferr(indb),tdphz(indb)],  ferrorBinEdges{e},  phzBins);
            end
        end        
    end
end


mmvjpdf =  zeros([1,numel(phzBinCenters),edc, 2]);
for y = 1:edc(1)
    for x = 1:edc(2)
        for z = 1:edc(3)
            for p = 1:numel(phzBinCenters),
                for e = 1:2,
                    try
                        [g, index] = unique(cumsum(mvjpdf(:,p,x,y,z,e)./sum(mvjpdf(:,p,x,y,z,e))),'last'); 
                        index(g==0|g==1) = [];
                        g(g==0|g==1) = [];
                        mmvjpdf(:,p,x,y,z,e) = interp1(g,ferrorBinCenters{e}(index),0.5);
                    end
                end
            end
        end
    end
end



for e = 1:2;
    for z = 1:3,
        figure();
        sax = reshape(tight_subplot(edc(1),edc(2),0.001,0,0)',[edc(2),edc(1)])';        
        for x = 1:edc(1),
            for y = 1:edc(2)
                axes(sax(x,y));                    
                imagesc(ferrorBinCenters{1},...
                        phzBinCenters(phzOrder)+2.*pi.*double(phzBinCenters(phzOrder)<0),...
                        imgaussfilt(mvjpdf(:,phzOrder,x,y,z,e),[2,0.2])');
                axis('xy');
                xlim([-300,300]); 
                set(gca(),'XTick',[])
                set(gca(),'YTick',[])
                Lines([],pi,'k');
                Lines(0,[],'k');
                colormap('jet');
            end
        end      
        set(gcf,'Units','centimeters');
        set(gcf,'Position',[(e-1)*7,z*3+10,7,3])
    end
end





%    'hba'    'hvl'    'hvf'    'hbp'
figure();
hold('on');
cmap = 'gcb';
offset = [-50:50:50];
for g = 1:3,
    subplot(3,1,4-g);    
    hold('on');
    for j = 1:varBinCnt,
        for i = 1:varBinCnt,
            plot(sq(mmvjpdf(1,phzOrder(fliplr(2:7)),j,i,g,2))+offset(j),...
                 sq(mmvjpdf(1,phzOrder(fliplr(2:7)),j,i,g,1)),[cmap(i),'+-'],'LineWidth',2)
            plot(sq(mmvjpdf(1,phzOrder(fliplr(5)),  j,i,g,2))+offset(j),...
                 sq(mmvjpdf(1,phzOrder(fliplr(5))  ,j,i,g,1)),[cmap(i),'o'],'LineWidth',2)
        end
    end
    xlim([offset(1)-diff(offset(1:2)),offset(end)+diff(offset(1:2))]);
    ylim([-25,125])
    Lines(-offset,[],'k');
    Lines([],0,'k');
    Lines(offset,[],'k');
    Lines(0,[],'k');
    set(gca,'XTick',offset);
    set(gca,'XTickLabel', vBinLbl{1});
    ylabel(vBinLbl{3}{g});
end
set(gcf,'Units','centimeters');
set(gcf,'Position',[(e)*7,3+10,10,9])

%%%>>>


figure,
hist2([-dfhrvf(ind&ismember(dtind,[3,4,5])),...
        dfhrvl(ind&ismember(dtind,[3,4,5]))],...
        linspace(-20,120,100),linspace(-120,120,100));
caxis([0,2000])