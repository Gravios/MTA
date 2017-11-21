version = '02';

% Figure Settings ----------------------------------------------------------------------
OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'MjgEd2016/manuscript/Figures/Suplementary';
mkdir(fullfile(OwnDir,FigDir));
% --------------------------------------------------------------------------------------


% LOAD Trial objects
% LOAD behavioral state collection
% LOAD marker position data
% PREFILTER xyz data
sessionList = get_session_list('ncp');
numTrials = numel(sessionList);
sessionList = mat2cell(sessionList,1,ones([1,numTrials]));
Trials = cf(@(t) MTATrial.validate(t)                    ,sessionList);
stc    = cf(@(t) t.load('stc','msnn_ppsvd')              ,Trials);
%xyz    = cf(@(t) preproc_xyz(t,'SPLINE_SPINE_HEAD_EQI')  ,Trials);
xyz    = cf(@(t) preproc_xyz(t)                          ,Trials);
         cf(@(x) x.filter('RectFilter')                  ,xyz);

% CREATE Virtual marker, on axis orthogonal to the transverse plane of the rat's head
nz  = cf(@(x)   -cross(x(:,'head_back',:)-x(:,'hcom',:),x(:,'head_left',:)-x(:,'hcom',:)), xyz);
nz  = cf(@(n)   bsxfun(@rdivide,n,sqrt(sum((n).^2,3))), nz);
nm  = cf(@(x,n) n.*20+x(:,'hcom',:)                   , xyz, nz);
      cf(@(x,n) x.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},n),  xyz,nm);
% CREATE Virtual marker, a low pass filtered copy of the head's center of mass
      cf(@(x)   x.addMarker('flhcom',[.7,1,.7],...
                            {{'head_back','head_front',[0,0,1]}},...
                            ButFilter(x(:,'hcom',:),3,2./(x.sampleRate/2),'low')),xyz);

% COMPUTE intermarker angles
% LOAD nasal cavity pressure sensor data
% LOAD rhythmic head motion 
% LOAD head pitch feature
% LOAD body refreneced position and motion
% LOAD head pitch and map to reference session
ang = cf(@(t,x) create(MTADang,t,x)    ,Trials,xyz);
ncp = cf(@(t,s) fet_ncp(t,s.ncpChannel),Trials,sessionList);
rhm = cf(@(t)   fet_rhm(t)             ,Trials);
rhp = cf(@(t)   fet_rhp(t)             ,Trials);
brf = cf(@(t)   fet_bref(t)            ,Trials);
pch = cf(@(t)   fet_HB_pitch(t)        ,Trials);
      cf(@(p,t) map_to_reference_session(p,t,'Ed05-20140529.ont.all'),pch,Trials);


% FILTER body referenced position and motion
% DETECT steps
fbrf = cf(@(b) b.copy(),  brf);
       cf(@(f) f.filter('ButFilter',3,[1.2,7],'bandpass'), fbrf);
steps = cf(@(f)   LocalMinima(-abs(f.data(:,17)),10,-3) ,fbrf);
steps = cf(@(s,p) s(WithinRanges(s,p{'w'}.data))  ,steps,stc);


% COMPUTE longitudinal rhythmic rotational motion of the head
rtm = cf(@(r)     r.copy()                                           , rhm);
%cf(@(r,a)   set(r,'data',circ_dist(a(:,'flhcom','htx',2),a(:,'flhcom','head_back',2))), rtm,ang);
cf(@(r,a)   set(r,'data',a(:,'flhcom','htx',3)), rtm,ang);
%cf(@(r,a)   set(r,'data',a(:,'spine_upper','hcom',3)), rtm,ang);
%cf(@(r,a)   set(r,'data',circ_dist(a(:,'flhcom','htx',2),a(:,'flhcom','head_front',2))), rtm,ang);
%cf(@(r,a,p) set(r,'data',unwrap(circ_dist(a.data,p.data))), rtm,prhm,prhp);
cf(@(r)     set(r,'data',[0;diff(RectFilter(diff(RectFilter(r.data))));0]), rtm);


% COMPUTE rhm power

[rhmp,fs,ts] = cf(@(t)   fet_rhm(t,[],'mtchglong'), Trials);
               cf(@(r,f) set(r,'data',mean(r(:,6<f&f<13),2,'omitnan')), rhmp,fs)
               cf(@(r,x) r.resample(x)            , rhmp,xyz);
for s = 1:numTrials, rhmp{s}.data(rhmp{s}.data<1e-15) = 1e-15; end
               cf(@(r)   set(r,'data',log10(r.data)), rhmp);

% COMPUTE phases
prhm = cf(@(r) r.phase([5,12]), rhm);
prhp = cf(@(r) r.phase([5,12]), rhp);
prtm = cf(@(r) r.phase([5,12]), rtm);
pncp = cf(@(n) n.phase([5,12]), ncp);

% CHECK reliability of signals
% $$$ figure();
% $$$ for s = 1:numTrials,
% $$$     subplot(numTrials,1,s);
% $$$     hold('on');
% $$$     plot(nunity(ncp{s}.data));
% $$$     plot(nunity(rhm{s}.data));
% $$$     plot(nunity(rhp{s}.data));
% $$$     ylim([-10,10]);
% $$$ end


[mins,vals] = cf(@(r)   LocalMinima(r.data(:,1),30,-0.01), rhm);
%[mins,vals] = cf(@(r)   LocalMinima(r.data(:,1),30,-1000), ncp);
mins = cf(@(m,s) m(WithinRanges(m,[s{'w+n+p+r'}.data])), mins,stc);

for s = 1:numTrials, mins{s}([1:10,end-10:end]) = []; end



hfig = figure();
% TOP ROW - ncp vs features
subplot(231); % pitch X phaseDiff(ncp,rhm)
outNcpRhm = cf(@(a,m,pn,pr)...
               histcounts2(circ_dist(pn(m),pr(m)),a(m,3),...
                           linspace(-pi,pi,50),linspace(-pi/2,pi/2,50)),...
               pch,mins,pncp,prhm);
imagesc(linspace(-pi,pi,50),linspace(-pi/2,pi/2,50),sum(cat(3,outNcpRhm{:}),3)');
axis('xy');title('pitch X phaseDiff(ncp,rhm)');

subplot(232); % pitch X phaseDiff(ncp,rhp)
outNcpRhp = cf(@(a,m,pn,pr)...
               histcounts2(circ_dist(pn(m),pr(m)),a(m,3),...
                           linspace(-pi,pi,50),linspace(-pi/2,pi/2,50)),...
               pch,mins,pncp,prhp);
imagesc(linspace(-pi,pi,50),linspace(-pi/2,pi/2,50),sum(cat(3,outNcpRhp{:}),3)');
axis('xy');title('pitch X phaseDiff(ncp,rhp)');

subplot(233); % ncp X rtm
outRhmRtm = cf(@(a,m,pn,pr)...
               histcounts2(circ_dist(pn(m),pr(m)),a(m,3),...
                           linspace(-pi,pi,50),linspace(-pi/2,pi/2,50)),...
               pch,mins,pncp,prtm);
imagesc(linspace(-pi,pi,50),linspace(-pi/2,pi/2,50),sum(cat(3,outRhmRtm{:}),3)');
axis('xy');title('pitch X phaseDiff(ncp,rtm)');


% BOTTOM ROW - pitch X phase(feature) - phase(feature) )
subplot(234); % pitch X phaseDiff(rhm,rhp)
outRhmRhp = cf(@(a,m,pn,pr)...
               histcounts2(circ_dist(pn(m),pr(m)),a(m,3),...
                           linspace(-pi,pi,50),linspace(-pi/2,pi/2,50)),...
               pch,mins,prhm,prhp);
imagesc(linspace(-pi,pi,50),linspace(-pi/2,pi/2,50),sum(cat(3,outRhmRhp{:}),3)');
axis('xy');title('pitch X phaseDiff(rhm,rhp)')

subplot(235); % pitch X phaseDiff(rhm,rhp)
outRhmRhp = cf(@(a,m,pn,pr)...
               histcounts2(circ_dist(pn(m),pr(m)),a(m,3),...
                           linspace(-pi,pi,50),linspace(-pi/2,pi/2,50)),...
               pch,mins,prtm,prhm);
imagesc(linspace(-pi,pi,50),linspace(-pi/2,pi/2,50),sum(cat(3,outRhmRhp{:}),3)');
axis('xy');title('pitch X phaseDiff(rtm,rhm)')

subplot(236); % pitch X phaseDiff(rtm,rhp)
outRhmRhp = cf(@(a,m,pn,pr)...
               histcounts2(circ_dist(pn(m),pr(m)),a(m,3),...
                           linspace(-pi,pi,50),linspace(-pi/2,pi/2,50)),...
               pch,mins,prtm,prhp);
imagesc(linspace(-pi,pi,50),linspace(-pi/2,pi/2,50),sum(cat(3,outRhmRhp{:}),3)');
axis('xy');title('pitch X phaseDiff(rtm,rhp)')

ForAllSubplots('caxis([0,100])')

FigName = ['pitch_phaseDiff_ncpXfet'];
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,version,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,version,'.png']));







% Tendency for zero phase offset between rhm and ncp for higher pitches
% Check phaseDif variance per cycle
dpnr = {};
for s = 1:numTrials,
    for m = 1:numel(mins{s}),
        dpnr{s}(m) = circ_var(circ_dist(prhm{s}(mins{s}(m):mins{s}(m)+18),prhp{s}(mins{s}(m):mins{s}(m)+18)));
    end
end        
figure();
outRhmRhp = cf(@(a,m,dpnr)...
               histcounts2(log10(dpnr)',a(m,3),...
                           linspace(-4,1,50),linspace(-pi/2,pi/2,50)),...
               pch,mins,dpnr);
imagesc(linspace(-4,1,50),linspace(-pi/2,pi/2,50),sum(cat(3,outRhmRhp{:}),3)');
axis('xy');title('pitch X phaseDiff(rhm,rhp)')






figure();
for i = 1:2:1500,

    plot(unwrap(circ_dist(prhm{1}(mins{1}(i):mins{1}(i+1)),prhp{1}(mins{1}(i):mins{1}(i+1)))))
    end
    pause(0.1)
end

figure,plot(unwrap(circ_dist(prhm{1}.data,prhp{1}.data)))



% phase{rhm X rhp} at steps
figure();
outNcpRhmRtm = cf(@(m,pn,pr)...
               histcounts2(pn(m),pr(m),linspace(-pi,pi,50),linspace(-pi,pi,50)),...
               steps,prhm,prhp);
imagesc(linspace(-pi,pi,50),linspace(-pi,pi,50),sum(cat(3,outNcpRhmRtm{:}),3)');
axis('xy');title('phase{rhm X rhp} at steps')

% phase{rhm X rtm} at steps
figure();
outNcpRhmRtm = cf(@(m,pn,pr)...
               histcounts2(pn(m),pr(m),linspace(-pi,pi,50),linspace(-pi,pi,50)),...
               steps,prhm,prtm);
imagesc(linspace(-pi,pi,50),linspace(-pi,pi,50),sum(cat(3,outNcpRhmRtm{:}),3)');
axis('xy');title('phase{rhm X rhp} at steps')

figure();
outNcpRhmRtm = cf(@(m,pn,pr)...
               histcounts2(pn(m),pr(m),linspace(-pi,pi,50),linspace(-pi,pi,50)),...
               mins,prhm,prhp);
imagesc(linspace(-pi,pi,50),linspace(-pi,pi,50),sum(cat(3,outNcpRhmRtm{:}),3)');
axis('xy');title('phase{rhm X rhp} at steps')


figure();
outNcpRhmRtm = cf(@(pn,pr)...
               histcounts2(pn(:),pr(:),linspace(-pi,pi,50),linspace(-pi,pi,50)),...
               prhm,prhp);
imagesc(linspace(-pi,pi,50),linspace(-pi,pi,50),sum(cat(3,outNcpRhmRtm{:}),3)');
axis('xy');title('phase{rhm X rhp} at steps')




figure();
outNcpRhmRtm = cf(@(pn,pr)...
               histcounts2(pn(:),pr(:),linspace(-pi,pi,50),linspace(-pi,pi,50)),...
               prhm,prhp);
imagesc(linspace(-pi,pi,50),linspace(-pi,pi,50),sum(cat(3,outNcpRhmRtm{:}),3)');
axis('xy');title('phase{rhm X rhp} at steps')


wrtm = nunity([rhm{1}.data,rhp{1}.data,ncp{1}.data]);
%wrtm(nniz(rtm)) = WhitenSignal(rtm(nniz(rtm)));
[ys,fs,ts] = mtchglong(wrtm,2^8,rtm{1}.sampleRate,2^7,2^7*0.875,[],[],[],[1,20]);

figure();
subplot(611);imagesc(ts,fs,log10(ys(:,:,1,1))');axis('xy');colormap('jet');caxis([-3,-1.5]);
subplot(612);imagesc(ts,fs,abs(ys(:,:,1,3))');       axis('xy');colormap('jet');caxis([0.5,1]);
subplot(613);imagesc(ts,fs,log10(ys(:,:,2,2))');axis('xy');colormap('jet');caxis([-2,-.5]);
subplot(614);imagesc(ts,fs,abs(ys(:,:,2,3))');       axis('xy');colormap('jet');caxis([0.5,1]);
subplot(615);imagesc(ts,fs,log10(ys(:,:,3,3))');axis('xy');colormap('jet');caxis([-3,0]);
subplot(616);plotSTC(stc{1},1);
linkaxes(findobj(gcf,'Type','Axes'),'x');


cmap = jet(100);

cind = cf(@(a,m) discretize(clip(a(m,3),-pi/2,pi/3),linspace(-pi/2,pi/3,100)), pch,mins);

cs = cf(@(c,m) m(c,:), cind,repmat({cmap},[1,numTrials]));

hfig = figure();
hold('on');
for s=1:numTrials,
    scatter(circ_dist(prhm{s}(mins{s}),prhp{s}(mins{s})),log10(abs(vals{s})),20,cs{s},'filled');
end
xlim([-pi,pi]);
ylim([-2,0]);
ylabel('log10(rhm peaks)');
xlabel('phaseDiff (rad)');
hcb = colorbar;
colormap(cmap);
caxis([-pi/2,pi/3]);
ylabel(hcb,'Pitch (rad)');
title('log10(rhm maxima) vs phaseDiff(rhm,rhp)');

FigName = ['phaseDiff_pitchXrhm_rhmPow'];
print(hfig,'-depsc2',fullfile(OwnDir,FigDir,[FigName,version,'.eps']));
print(hfig,'-dpng',  fullfile(OwnDir,FigDir,[FigName,version,'.png']));







% $$$ Trial = MTATrial.validate('Ed05-20140529.ont.all');
% $$$ % LOAD behavioral state collection
% $$$ stc = Trial.load('stc','msnn_ppsvd');
% $$$ % LOAD marker position data
% $$$ xyz = preproc_xyz(Trial);
% $$$ % PREFILTER xyz data
% $$$ xyz.filter('RectFilter');
% $$$ % CREATE Virtual marker, on axis orthogonal to the transverse plane of the rat's head
% $$$ nz  = -cross(xyz(:,'head_back',:)-xyz(:,'hcom',:),xyz(:,'head_left',:)-xyz(:,'hcom',:));
% $$$ nz  = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
% $$$ nm  = nz.*20+xyz(:,'hcom',:);
% $$$ xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
% $$$ % CREATE Virtual marker, a low pass filtered copy of the head's center of mass
% $$$ xyz.addMarker('flhcom',[.7,1,.7],...
% $$$               {{'head_back','head_front',[0,0,1]}},...
% $$$               ButFilter(xyz(:,'hcom',:),3,2./(xyz.sampleRate/2),'low'));
% $$$ % COMPUTE intermarker angles
% $$$ ang = create(MTADang,Trial,xyz);
% $$$ % LOAD nasal cavity pressure sensor data
% $$$ ncp = fet_ncp(Trial);
% $$$ % LOAD rhythmic head motion 
% $$$ rhm = fet_rhm(Trial);
% $$$ % COMPUTE longitudinal rhythmic rotational motion of the head
% $$$ rtm = rhm.copy();
% $$$ rtm.data = circ_dist(ang(:,'flhcom','htx',2),ang(:,'flhcom','head_front',2));
% $$$ rtm.filter('RectFilter');
% $$$ rtm.data = [0;diff(RectFilter(diff(rtm.data)));0];
% $$$ 
% $$$ % COMPUTE phases
% $$$ prhm = rhm.phase([5,12]);
% $$$ prtm = rtm.phase([5,12]);
% $$$ pncp = ncp.phase([5,12]);
% $$$ 
% $$$ 
% $$$ 
% $$$ wrtm = nunity([rhm.data,rtm.data,ncp.data]);
% $$$ %wrtm(nniz(rtm)) = WhitenSignal(rtm(nniz(rtm)));
% $$$ [ys,fs,ts] = mtchglong(wrtm,2^8,rtm.sampleRate,2^7,2^7*0.875,[],[],[],[1,20]);
% $$$ 
% $$$ figure();
% $$$ subplot(611);imagesc(ts,fs,log10(ys(:,:,1,1))');axis('xy');colormap('jet');caxis([-4,-0.5]);
% $$$ subplot(612);imagesc(ts,fs,abs(ys(:,:,1,3))');       axis('xy');colormap('jet');caxis([0.25,1]);
% $$$ subplot(613);imagesc(ts,fs,log10(ys(:,:,2,2))');axis('xy');colormap('jet');caxis([-3,0]);
% $$$ subplot(614);imagesc(ts,fs,abs(ys(:,:,2,3))');       axis('xy');colormap('jet');caxis([0.25,1]);
% $$$ subplot(615);imagesc(ts,fs,log10(ys(:,:,3,3))');axis('xy');colormap('jet');caxis([-3,0]);
% $$$ subplot(616);plotSTC(stc,1);
% $$$ linkaxes(findobj(gcf,'Type','Axes'),'x');
% $$$ 
% $$$ % LOAD bref
% $$$ bref = fet_bref(Trial);
% $$$ fbref= bref.copy();
% $$$ fbref.filter('ButFilter',3,[1.2,7],'bandpass');
% $$$ % DETECT steps
% $$$ mins = LocalMinima(fbref.data(:,17),10,-3);
% $$$ mins = mins(WithinRanges(mins,stc{'w'}.data));
% $$$ 
% $$$ pfbref = fbref.phase([5,12]);
% $$$ 
% $$$ figure,
% $$$ i = -18;
% $$$ plot(bref(mins-i,16),prhm(mins-i),'.')
% $$$ 
% $$$ plot(bref(mins,27),prhm(mins),'.')
% $$$ 
% $$$ figure,
% $$$ for i = 1:30,
% $$$     subplot(5,6,i);
% $$$     rose(circ_dist(pfbref(mins-i,19),prhm(mins)))
% $$$ end
% $$$ 
% $$$ 
% $$$ plot(bref(mins,16),circ_dist(pfbref(mins-i,26),prhm(mins)),'.')
% $$$ 
% $$$ figure,rose(circ_dist(prtm(mins,18),prhm(mins)))
% $$$ 
% $$$ figure,hist(circ_dist(pfbref(mins,27),prhm(mins)),30)
% $$$ 
% $$$ 
% $$$ 
% $$$ mins = LocalMinima(rhm.data(:,1),10,-0.1);
% $$$ mins = mins(WithinRanges(mins,[stc{'w+n+p+r'}.data]));
% $$$ 
% $$$ figure,plot(ang(mins,'head_back','head_front',2),circ_dist(prtm(mins),prhm(mins)),'.')
% $$$ 
% $$$ figure();
% $$$ histogram2(ang(mins,'head_back','head_front',2),circ_dist(pncp(mins),prhm(mins)),[50,50],'DisplayStyle','tile')
% $$$ histogram2(ang(mins,'head_back','head_front',2),circ_dist(prtm(mins),pncp(mins)),[50,50],'DisplayStyle','tile')
% $$$ histogram2(ang(mins,'head_back','head_front',2),circ_dist(prhm(mins),prtm(mins)),[50,50],'DisplayStyle','tile')
% $$$ 
% $$$ 
