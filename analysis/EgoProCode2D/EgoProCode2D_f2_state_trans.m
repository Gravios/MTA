


% ---------------------------------------------------------------------------
% Characterization of the relationship between the rat's radial position on
% the maze and the orientation of the rat relative to the center of the maze.

% SETUP 
opts.title.fontsize = {'FontSize',14};
opts.labels.fontsize = {'FontSize',14};

binz.hba.label = 'hba';
binz.hba.name = 'Head Body Angle'
binz.hba.edges = linspace(-1.5,1.5,20);
binz.hba.centers = mean([binz.hba.edges(1:end-1);binz.hba.edges(2:end)]);
binz.hba.count = numel(binz.hba.centers);

binz.hma.label = 'hma';
binz.hma.name = 'Head Maze Angle'
binz.hma.edges = linspace(-pi,pi,30);
binz.hma.centers = mean([binz.hma.edges(1:end-1);binz.hma.edges(2:end)]);
binz.hma.count = numel(binz.hma.centers);

binz.dst.label = 'dst';
binz.dst.name = 'Head Maze Distance';
binz.dst.edges = [0,200,300,500];
binz.dst.centers = mean([binz.dst.edges(1:end-1);binz.dst.edges(2:end)]);
binz.dst.count = numel(binz.dst.centers);

samplerate = dca{1}.sampleRate;

% COMPUTATION  
hout = zeros([binz.hma.count,binz.hba.count,binz.dst.count,numel(dca)]);
for t = 1:numel(dca);
    stc = Trials{tind(t)}.load('stc','msnn_ppsvd');
    % anteroposterior component of egocentric decoding
    fwd = nan([size(xyz{tind(t)},1),1]);  fwd(dca{t}.ind) = dca{t}.esax(:,1);
    % lateral component of egocentric decoding
    lat = nan([size(xyz{tind(t)},1),1]);  lat(dca{t}.ind) = dca{t}.esax(:,2);
    % theta phase
    phz = nan([size(xyz{tind(t)},1),1]);  phz(dca{t}.ind) = dca{t}.phz(:,1);
    % distance to maze center
    dst = nan([size(xyz{tind(t)},1),1]);  dst(dca{t}.ind) = dca{t}.hdist(:,1);
    % head body angle
    hba = nan([size(xyz{tind(t)},1),1]);  hba(dca{t}.ind) = dca{t}.hbang(:,1);
    % head angular velocity
    hva = nan([size(xyz{tind(t)},1),1]);  hva(dca{t}.ind) = dca{t}.hvang(:,1);
    stsT = nan([size(xyz{tind(t)},1),1]);
    stsT(dca{t}.ind) = double(dca{t}.stcm(:,1)==1);
    sss = double(dca{t}.stcm(:,1)==1);
    sss(sss==0) = nan;
    stsT(~isnan(stsT)) = sss;
    hma = nan([size(xyz{tind(t)},1),1]);
    headAngle = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
    headAngle = atan2(headAngle(:,2),headAngle(:,1));
    mazeAngle = sq(dca{t}.xyz(:,'hcom',[1,2]));
    mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
    hma(dca{t}.ind) = circ_dist( headAngle, mazeAngle);
    
    for id = 1:binz.dst.count
        ind =  binz.dst.edges(:,id) < dst&dst < binz.dst.edges(:,id+1);
        hout(:,:,id,t) = ...
            hist2([nonzeros(hma(ind).*stsT(ind)), ...
                   nonzeros(hba(ind).*stsT(ind))],...
                  binz.hma.edges,...
                  binz.hba.edges);
    end
end

% PLOTTING 
figure();
for id = 1:binz.dst.count
    subplot(binz.dst.count,1,id);
    dout = sum(hout(:,:,id,:),4);
    dout = bsxfun(@rdivide,dout,sum(dout));
    imagesc(binz.hma.centers,...
            binz.hba.centers,...
            dout');
    
    title({[binz.hma.label,' VS ',binz.hba.label],...
           num2str(binz.dst.edges(id:id+1),'proximity to center ≬[%d,%d]cm')},...
           opts.title.fontsize{:});
    xlabel([binz.hma.label, ' (rad)'],opts.labels.fontsize{:});
    ylabel([binz.hba.label, ' (rad)'],opts.labels.fontsize{:});    
    daspect([1,1,1]);
    axis('xy');
    colormap('jet')
    caxis([0,0.2])
    cax = colorbar()
    ylabel(cax,'Probability',opts.labels.fontsize{:});
end
set(gcf(),'PaperOrientation','portrait');
set(gcf(),'PaperType','A4');
% ---------------------------------------------------------------------------











t = 7;
samplerate = dca{1}.sampleRate;

stc = Trials{tind(t)}.load('stc','msnn_ppsvd');


%tp = Trials{tind(t)}.stc.get_state_transitions(Trials{tind(t)},{'walk','pause'},0.1,xyz{tind(t)});
tp = stc.get_state_transitions(Trials{tind(t)},{'pause','turn'},0.1,xyz{tind(t)});
tp = cat(1,tp,stc.get_state_transitions(Trials{tind(t)},{'walk','turn'},0.1,xyz{tind(t)}));

tp = stc.get_state_transitions(Trials{tind(t)},{'turn','walk'},0.1,xyz{tind(t)});


% anteroposterior component of egocentric decoding
fwd = nan([size(xyz{tind(t)},1),1]);
fwd(dca{t}.ind) = dca{t}.esax(:,1);
% lateral component of egocentric decoding
lat = nan([size(xyz{tind(t)},1),1]);
lat(dca{t}.ind) = dca{t}.esax(:,2);
% theta phase
phz = nan([size(xyz{tind(t)},1),1]);
phz(dca{t}.ind) = dca{t}.phz(:,1);
% distance to maze center
dst = nan([size(xyz{tind(t)},1),1]);
dst(dca{t}.ind) = dca{t}.hdist(:,1);
% head body angle
hba = nan([size(xyz{tind(t)},1),1]);
hba(dca{t}.ind) = dca{t}.hbang(:,1);
% head angular velocity
hva = nan([size(xyz{tind(t)},1),1]);
hva(dca{t}.ind) = dca{t}.hvang(:,1);

stsT = nan([size(xyz{tind(t)},1),1]);
stsT(dca{t}.ind) = double(dca{t}.stcm(:,1)==1);
sss = double(dca{t}.stcm(:,1)==1);
sss(sss==0) = nan;
stsT(~isnan(stsT)) = sss;





% sturcture of the search
%  - ¿ does the current hba side influence the direction of the turn ?
%      - account for location based geometric restrictions on turn (direction / probability)
%          - (object / wall) valence may out rank (maze center / distant cues) 
%  - ¿ What is the likelihood of the an isplateral turn given both hba and previous-turn sign ?

hvfl = fet_hvfl(Trials{tind(t)},samplerate);
fhvl = hvfl.copy();
fhvl.data = fhvl.data(:,2);
fhvl.filter('ButFilter',4,2,'low');

figure,plot(hvfl(:,2),'k'),hold('on');plot(fhvl.data,'r');


figure();
subplot(211);
plot([1:size(fhvl,1)]./samplerate,abs(fhvl.data),'r');
subplot(212);
plotSTC(stc,1);
linkx();

tp = ThreshCross(abs(fhvl.data),10,10);
tp = tp(1:3:end,:);


fhba = fet_hba(Trials{t}, samplerate);
fhba.filter('ButFilter',4,2,'low');

tx = nonzeros(ThreshCross(fhba.data,0,10)');

figure,histogram(log10(diff(tx)./samplerate),linspace(-2,2,100))
figure,histogram(nonzeros(tx-tx')./samplerate,linspace(-2,2,100))


figure();hist2(log10(abs(diff([tx,circshift(tx,-1)]./samplerate))),linspace(-2,2,40),linspace(-2,2,40),'')

txda = nonzeros((tx-tx').*double(~eye(numel(tx))));
txds = nonzeros((diff([0;tx])+diff([tx;size(fhba,1)]))'.*double(~eye(numel(tx))));


% the threshold gets timepoints where the head-body orientation changes between left and right.
% the difference between neighboring time points provides the durations in left or right states.




figure();
subplot(211);hist2([txda./samplerate,log10(abs(txds/samplerate))],linspace(-10,10,80),linspace(-1,2,40),'yprob')
caxis([0,0.1])
subplot(212);hist2([txda./samplerate,circshift(log10(abs(txds/samplerate)),1)],linspace(-10,10,80),linspace(-1,2,40),'yprob')
caxis([0,0.1])

% head maze angle [0: toward maze center, pi: toward maze wall]
hma = nan([size(xyz{tind(t)},1),1]);
headAngle = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
headAngle = atan2(headAngle(:,2),headAngle(:,1));
mazeAngle = sq(dca{t}.xyz(:,'hcom',[1,2]));
mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
hma(dca{t}.ind) = circ_dist( headAngle, mazeAngle);

mfun = @(beta,x) beta(2).*(sum(x.^2,2)).*sin(atan2(x(:,2),x(:,1)))+beta(1);
[Xea,Yea] = pol2cart(hma,dst);
clt = lat-mfun(beta,[Xea(:),Yea(:)]);


figure();
plot(lat.*stsT);
hold('on');
plot(clt.*stsT);
plot(hma.*200.*stsT);
plot(hba.*200.*stsT);

figure,plot(hma.*stsT,hba.*stsT,'.');




figure();
hold('on');
plot(dst,'k');                           % DISTANCE 
plot(lat-mfun(beta,[Xea(:),Yea(:)]),'c');% CORR LAT
plot(lat,'m');                           % NORM LAT
plot(mfun(beta,[Xea(:),Yea(:)]),'r');    % CORR FUN
hold('on'),plot(hba*50)                  % HBA


% Maze Border -> Circle
% Maze Center -> Small circle
% Subject -> Subject
% Line From MazeCenter to SubjectHead
% Line From SubjectNeck to SubjectNose
% Line From SubjectTail to SubjectNeck
% Line From SubjectNeck extexending previous




 
hwin =250;
p2l.tpr = round(mean(tp,2));
p2l.fwd = GetSegs(fwd, p2l.tpr-hwin,2*hwin);
p2l.lat = GetSegs(lat, p2l.tpr-hwin,2*hwin);
p2l.clt = GetSegs(clt, p2l.tpr-hwin,2*hwin);
p2l.phz = GetSegs(phz, p2l.tpr-hwin,2*hwin);
p2l.dst = GetSegs(dst, p2l.tpr-hwin,2*hwin);
p2l.hba = GetSegs(hba, p2l.tpr-hwin,2*hwin);
p2l.hva = GetSegs(hva, p2l.tpr-hwin,2*hwin);
p2l.hma = GetSegs(hma, p2l.tpr-hwin,2*hwin);


mdist = 200;
rind =   dst(p2l.tpr)>mdist;
%         | ~WithinRanges(-hma(p2l.tpr),[pi/4,pi*3/4]);
%          | ~WithinRanges(abs(hma(p2l.tpr)),[0,pi*1/4]);
%         | ~WithinRanges(hma(p2l.tpr),[pi/4,pi*3/4]);
%          | ~WithinRanges(abs(hma(p2l.tpr)),[0,pi/4]);

p2l.fwd(:,rind) = [];
p2l.lat(:,rind) = [];
p2l.clt(:,rind) = [];
p2l.phz(:,rind) = [];
p2l.hba(:,rind) = [];
p2l.hva(:,rind) = [];
p2l.hma(:,rind) = [];
p2l.tpr(rind) = [];

p2l.fwd(~WithinRanges(p2l.phz,[4.5,5.5])) = nan;
p2l.lat(~WithinRanges(p2l.phz,[4.5,5.5])) = nan;
p2l.clt(~WithinRanges(p2l.phz,[4.5,5.5])) = nan;

%out = sq(mean(reshape(p2l.fwd,[50,20,size(p2l.lat,2)]),2,'omitnan'));
out = sq(mean(reshape(p2l.clt,[2*hwin/20,20,size(p2l.clt,2)]),2,'omitnan'));

%bind = round(hwin*2.*[0.45,0.5]);
bind = round(hwin*2.*[0.5,0.6]);
bind = bind(1):bind(2);
[hbs,si] = sort(mean(p2l.hma(bind,:),'omitnan'));
ots = ([1:2*hwin/20]-12.5)/12.5;
bts = linspace(-hwin/xyz{tind(t)}.sampleRate,...
               hwin/xyz{tind(t)}.sampleRate,...
               2*hwin);

ny = 7;
p = 1;
figure();
subplot(ny,1,p);p=p+1;
    plot(1:size(p2l.clt,2),hbs);
    xlim([1,size(p2l.clt,2)]);
    Lines([],0,'k');
subplot(ny,1,p);p=p+1;
    imagescnan({1:numel(si),ots,out(:,si)},'colorLimits',[-200,200],'colorMap',@jet);
subplot(ny,1,p);p=p+1; % Mean lateral 
    hold('on');
    plot(1:size(out,2),mean(out(1:10,si)-18,'omitnan'));
    plot(1:size(out,2),mean(out(16:25,si)-18,'omitnan'));
    %plot(1:size(out,2),mean(out(16:25,si),'omitnan')-mean(out(1:10,si),'omitnan'),'m');
    ylim([-200,200]);
    Lines([],0,'k');
subplot(ny,1,p);p=p+1;
    fh = find(hbs<0);
    sh = find(hbs>0);
    hold('on');
    plot(ots,mean(out(:,si(fh)),2,'omitnan')-18);
    plot(ots,mean(out(:,si(sh)),2,'omitnan')-18);
    ylim([-150,150]);
    Lines([],0,'k');
subplot(ny,1,p);p=p+1;
    hold('on');
    plot(bts,mean(p2l.hba(:,si(fh)),2,'omitnan'),'b');
    plot(bts,mean(p2l.hba(:,si(fh)),2,'omitnan')+std(p2l.hba(:,si(fh)),[],2,'omitnan'),'m');
    plot(bts,mean(p2l.hba(:,si(fh)),2,'omitnan')-std(p2l.hba(:,si(fh)),[],2,'omitnan'),'m');
    plot(bts,mean(p2l.hba(:,si(sh)),2,'omitnan'),'r');
    plot(bts,mean(p2l.hba(:,si(sh)),2,'omitnan')-std(p2l.hba(:,si(sh)),[],2,'omitnan'),'g');    
    plot(bts,mean(p2l.hba(:,si(sh)),2,'omitnan')+std(p2l.hba(:,si(sh)),[],2,'omitnan'),'g');    
    Lines([],0,'k');
subplot(ny,1,p);p=p+1;
    hold('on');
    plot(bts,mean(p2l.hva(:,si(fh)),2,'omitnan'));
    plot(bts,mean(p2l.hva(:,si(sh)),2,'omitnan'));
    Lines([],0,'k');
subplot(ny,1,p);p=p+1;
    hold('on');
    plot(bts,mean(p2l.hma(:,si(fh)),2,'omitnan'));
    plot(bts,mean(p2l.hma(:,si(sh)),2,'omitnan'));
    ylim([-pi,pi])
    Lines([],0,'k');



figure()
    hold(gca(),'on');
    histogram(mean(out(16:25,si(fh)),'omitnan')-mean(out(1:10,si(fh)),'omitnan'),linspace(-150,150,20),'FaceColor','r','FaceAlpha',0.3);
    histogram(mean(out(16:25,si(sh)),'omitnan')-mean(out(1:10,si(sh)),'omitnan'),linspace(-150,150,20),'FaceColor','g','FaceAlpha',0.3);


[h,p,c,s] = ttest2(mean(out(16:25,si(fh)),'omitnan')-mean(out(1:10,si(fh)),'omitnan'),mean(out(16:25,si(sh)),'omitnan')-mean(out(1:10,si(sh)),'omitnan'))



figure,plot(p2l.hma(:),p2l.clt(:),'.')



dind =   abs(p2l.hba(:))<1.2 ...
       & nniz(p2l.hba(:))    ...
       & nniz(p2l.clt(:));

figure,plot(p2l.hma(:),p2l.lat(:),'.')

figure();
subplot(221);
scatter(p2l.hba(dind),p2l.clt(dind),10,p2l.hva(dind),'filled');
colormap('jet');
caxis([-.5,.5]);

subplot(222);
scatter(p2l.hva(dind),p2l.clt(dind),10,p2l.hba(dind),'filled')
colormap('jet');
xlim([-.5,.5]);

subplot(223);
scatter(p2l.hba(dind),p2l.hva(dind),10,p2l.clt(dind),'filled')
colormap('jet');
ylim([-.5,.5]);

subplot(224);
scatter(p2l.hba(dind),p2l.hma(dind),10,p2l.clt(dind),'filled')
colormap('jet');
ylim([-pi,pi]);

[rho,pv ] = corr(p2l.hba(dind), p2l.clt(dind))


[rho,pv ] = corr(p2l.hva(dind), p2l.clt(dind))

[rho,pv ] = corr(p2l.hma(dind), p2l.clt(dind))

[rho,pv ] = corr(p2l.hma(dind), p2l.hba(dind))


figure();
scatter(p2l.hma(dind),p2l.clt(dind),10,p2l.hba(dind),'filled');
colormap('jet');
caxis([-1,1]);



figure,hist2([p2l.hba(dind),p2l.hma(dind)],linspace(-1.2,1.2,20),linspace(-pi,pi,30));


figure,hist2([p2l.hba(dind)],linspace(-1.2,1.2,20),

ihma = discretize(p2l.hma(dind),linspace(-pi,pi,30));
ihba = discretize(p2l.hba(dind),linspace(-1.2,1.2,20));

mout = accumarray([ihba,ihma],p2l.clt(dind),[20,30],@mean);


figure();
imagesc(mout');
colormap('jet');
axis('xy');
caxis([-200,200])









ahma = [];
aclt = [];
afwd = [];
alat = [];
ahba = [];
ahva = [];
aphz = [];
adst = [];



samplerate = dca{1}.sampleRate;
mfun = @(beta,x) beta(2).*(sum(x.^2,2)).*sin(atan2(x(:,2),x(:,1)))+beta(1);

for t = 1:numel(dca);
    stc = Trials{tind(t)}.load('stc','msnn_ppsvd');
    % anteroposterior component of egocentric decoding
    fwd = nan([size(xyz{tind(t)},1),1]);   fwd(dca{t}.ind) = dca{t}.esax(:,1);
    % lateral component of egocentric decoding
    lat = nan([size(xyz{tind(t)},1),1]);   lat(dca{t}.ind) = dca{t}.esax(:,2);
    % theta phase
    phz = nan([size(xyz{tind(t)},1),1]);   phz(dca{t}.ind) = dca{t}.phz(:,1);
    % distance to maze center
    dst = nan([size(xyz{tind(t)},1),1]);   dst(dca{t}.ind) = dca{t}.hdist(:,1);
    % head body angle
    hba = nan([size(xyz{tind(t)},1),1]);   hba(dca{t}.ind) = dca{t}.hbang(:,1);
    % head angular velocity
    hva = nan([size(xyz{tind(t)},1),1]);   hva(dca{t}.ind) = dca{t}.hvang(:,1);
    % head maze angle [0: toward maze center, pi: toward maze wall]
    hma = nan([size(xyz{tind(t)},1),1]);
    headAngle = sq(dca{t}.xyz(:,'nose',[1,2])-dca{t}.xyz(:,'hcom',[1,2]));
    headAngle = atan2(headAngle(:,2),headAngle(:,1));
    mazeAngle = sq(dca{t}.xyz(:,'hcom',[1,2]));
    mazeAngle = atan2(mazeAngle(:,2),mazeAngle(:,1));
    hma(dca{t}.ind) = circ_dist( headAngle, mazeAngle);
    % CORRECTED LAT    
    [Xea,Yea] = pol2cart(hma,dst);
    clt = lat-mfun(beta,[Xea(:),Yea(:)]);
    % THETA STATE
    stsT = nan([size(xyz{tind(t)},1),1]);
    stsT(dca{t}.ind) = double(dca{t}.stcm(:,1)==1 ...
                              & (dca{t}.stcm(:,3)==3 ...
                                 |dca{t}.stcm(:,4)==4 ...
                                 |dca{t}.stcm(:,5)==5));
    sss = double(dca{t}.stcm(:,1)==1);
    sss(sss==0) = nan;
    stsT(~isnan(stsT)) = sss;
    % CONCAT 
    ahma = cat(1,ahma,nonzeros(hma.*stsT));
    aclt = cat(1,aclt,nonzeros(clt.*stsT));
    afwd = cat(1,afwd,nonzeros(fwd.*stsT));
    alat = cat(1,alat,nonzeros(lat.*stsT));
    ahba = cat(1,ahba,nonzeros(hba.*stsT));
    ahva = cat(1,ahva,nonzeros(hva.*stsT));    
    aphz = cat(1,aphz,nonzeros(phz.*stsT));
    adst = cat(1,adst,nonzeros(dst.*stsT));
end





figure();
dbound = [0,150,300,450];
ny = numel(dbound)-1;

for d = 1:ny
    dind =   abs(ahba(:))<1.2 ...
             & nniz(ahba(:))    ...
             & nniz(aclt(:)) ...
             & adst<dbound(d+1) & adst>dbound(d) ...
             & WithinRanges( aphz, [4.5,5.5]);

    ihma = discretize(ahma(dind),linspace(-pi,pi,12));
    ihba = discretize(ahba(dind),linspace(-1.2,1.2,8));


    vclt = alat(dind);
    nind = nniz(ihba) & nniz(ihma) & nniz(vclt);
    mout = accumarray([ihba(nind),ihma(nind)],vclt(nind),[7,11],@mean);

    subplot2(ny,2,d,1);
        cout = hist2([ahba(dind),ahma(dind)],...
              linspace(-1.2,1.2,8),...
              linspace(-pi,pi,12),...
                     '');
        imagesc(linspace(-1.2,1.2,8),linspace(-pi,pi,12),cout');
        axis('xy');
        set(gca(),'YTick',[-pi,-pi/2,0,pi/2,pi]);
        set(gca(),'YTickLabels',{'Inbound','CW','Outbound','CCW'});
    
    subplot2(ny,2,d,2);

        mout(cout<50) = nan;
        imagescnan({linspace(-1.2,1.2,7),...
                    linspace(-pi,pi,11),...
                    (mout-18)'},...
                   'colorLimits',[-150,150],   ...
                   'colorMap',@jet);        
        axis('xy');
        
        xlabel([binz.hba.name,' (rad)'])
        ylabel([binz.hma.name,' (rad)'])
        title({[binz.hba.label,' vs ', binz.hma.label],...
               num2str(dbound(d:d+1),'Maze to Rat Dist [%d, %d]')});
end 

d = 3;
STATS =[];
dind =   abs(ahba(:))<1.2 ...
         & nniz(ahba(:))    ...
         & nniz(aclt(:)) ...
         & adst<dbound(d+1) & adst>dbound(d) ...
         & WithinRanges( aphz, [4.5,5.5]) ...
         & randn(size(ahba))>0;


[B,BINT,R,RINT,STATS(end+1,:)] = regress(aclt(dind),[ones([sum(dind),1]),ahba(dind),ahma(dind)]);
[B,BINT,R,RINT,STATS(end+1,:)] = regress(aclt(dind),[ones([sum(dind),1]),ahma(dind)]);
[B,BINT,R,RINT,STATS(end+1,:)] = regress(R,[ones([sum(dind),1]),ahba(dind)]);
STATS

%    1.6% ofthe variance


[B,BINT,R,RINT,STATS] = regress(aclt(dind),[ones([sum(dind),1]),ahba(dind),ahma(dind),adst(dind)]);
STATS

[B,BINT,R,RINT,STATS] = regress(aclt(dind),[ones([sum(dind),1]),ahba(dind),adst(dind)]);
STATS

[B,BINT,R,RINT,STATS] = regress(ahba(dind),[ones([sum(dind),1]),ahma(dind)]);
STATS


[B,BINT,R,RINT,STATS] = regress(aclt(dind),[ones([sum(dind),1]),ahba(dind)]);
STATS

[B,BINT,R,RINT,STATS] = regress(aclt(dind),[ones([sum(dind),1]),ahma(dind)]);
STATS

