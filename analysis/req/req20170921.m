% req20170921 ---------------------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: lateral motions of head and body synchrony during locomotion
%  Bugs: NA


Trial = MTATrial.validate('Ed05-20140529.ont.all');
stc = Trial.load('stc','hand_labeled_rev1_Ed');
xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQI');
features = fet_bref(Trial);
features.data = features.data(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30]);
fetHB = fet_HB_angvel(Trial);



fxyz = xyz.copy();
fxyz.filter('ButFilter',3,2.4,'low');



hbfet = features.copy();
hbfet.data = unity([hbfet.data(:,[21]),fetHB(:,2)]);
hbfet.data = hbfet.data(:,[21,25]);



specParm =    struct('nFFT',2^7,'Fs',features.sampleRate,...
                           'WinLength',2^6,'nOverlap',2^6*.875,...
                           'FreqRange',[1,10]);


[ys,fs,ts] = fet_spec(Trial,hbfet,'mtcsdglong',false,'defspec',specParm);

wper = stc{'w',ys.sampleRate};
vxy = fxyz.vel([1,2],[1,2]);
vxy.resample(ys);
vxy.data(vxy.data<1e-3) = nan;
%vxy.data = log10(vxy.data);


chr = [];
mvx = [];
mys = [];
man = [];
for w = wper.data',
    chr(end+1,:) = nanmean(abs(ys(w,:,1,2)))./nanmean(sqrt(ys(w,:,1,1).*ys(w,:,2,2)));
    mvx(end+1,:) = mean(vxy(w,:));
    mys(end+1,:) = nanmean(ys(w,:,2,2));
    %man(end+1,:) = circ_mean(fetHBp(w,:));
end


% FILTER by walk duration
wdur = diff(wper.data,1,2);
ind = wdur>=14;

out = [];
for f = 1:numel(fs),
    out(:,f) = histc(chr(ind,f),linspace(0,1,30));
end
figure,imagesc(1:30,fs,out');


figure,hist(chr(ind,5),linspace(0,1,40));
figure,plot(chr(ind,4),mvx(ind,2),'.');
figure,plot(chr(ind,5),log10(mys(ind,15))./log10(mys(ind,5)),'.');
figure,plot(chr(ind,5),man(ind,2),'.');


%% Multi Session 
OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_4/parts';
mkdir(fullfile(OwnDir,FigDir));


Trials   = af(@(t)  MTATrial.validate(t),  get_session_list('hand_labeled'));
stc      = cf(@(t)  t.load('stc'),         Trials);
stca     = cf(@(t)  t.load('stc',[t.stc.mode,'_SVDTRAJADJ']),         Trials);
xyz      = cf(@(t)  preproc_xyz(t,'SPLINE_SPINE_HEAD_EQI'),  Trials);

           
features = cf(@(t)  fet_bref(t),         Trials);
           cf(@(f)  set(f,'data',f(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30])),features);

lfet   = cf(@(f)   f.copy(),                           features);
         cf(@(f,p)   set(f,'data',[f.data(:,21:25)]),    lfet,features);
           
%[fetp]    = cf(@(t)  fet_HB_angvel(t),         Trials);           

for featureSet = {'fet_HB_pitchvel','fet_HB_pitch','fet_HB_angvel'},
    [hbfet,fett,fetd]    = cf(@(t,f)  feval(f,t),         Trials,repmat(featureSet(1),[1,numel(Trials)]));
    %[hbfet,fett,fetd]    = cf(@(t)  fet_HB_pitch(t),         Trials);           
    %[hbfet,fett,fetd]    = cf(@(t)  fet_HB_pitchvel(t),         Trials);
    %[hbfet,fett,fetd]    = cf(@(t)  fet_HB_angvel(t),         Trials);
% $$$ hbfet   = cf(@(f)   f.copy(),                           features);
% $$$           cf(@(f,p)   set(f,'data',[f.data(:,25),p.data(:,3)]),    hbfet,fetp);
% $$$ hbfet{1}.label = 'hbref_azvel';
%          cf(@(f)   set(f,'data',f.data(:,[21,25])),    hbfet);
%          cf(@(f)   set(f,'data',[f.data(:,[21]),flipud(f.data(:,[21,25]))]),    hbfet);

    
    specParm = repmat({struct('nFFT',2^9,'Fs',features{1}.sampleRate,...
                              'WinLength',2^8,'nOverlap',2^8*0.875,...
                              'FreqRange',[1,20])},...
                      [1,numel(Trials)]);
    %[ys,fs,ts] = cf(@(t,f,p)  fet_spec(t,f,'mtfft',false,'defspec',p),  Trials,hbfet,specParm);
    [ys,fs,ts] = cf(@(t,f,p)  fet_spec(t,f,'mtcsdglong',false,'defspec',p),  Trials,hbfet,specParm);

% $$$ fetp     = cf(@(t)   fet_HB_pitch(t),                    Trials);
% $$$            cf(@(f,t) f.map_to_reference_session(t,'jg05-20120317.cof.all'), fetp, Trials);
% $$$            cf(@(f,y) f.resample(y),                      fetp, ys);           
    
    fxyz     = cf(@(x)   x.copy(),                           xyz     );
    cf(@(x)   x.filter('ButFilter',3,2.4,'low'),  fxyz    );
    cf(@(x,y) x.resample(y),                      fxyz, ys);

    vxy      = cf(@(x)   x.vel([1,2],[1,2]),   fxyz);
    for s = 1:numel(vxy), vxy{s}.data(vxy{s}.data<1e-3) = 1e-3; end
    ang      = cf(@(t,x)  create(MTADang,t,x),  Trials, fxyz);


    for sts = 'wnr',
        wper     = cf(@(s,y) s{sts,y.sampleRate},  stc,  ys);
        %wper     = cf(@(s,y) s{'r',y.sampleRate},  stca,  ys);

        for fetInds = [1,1,2;2,3,3],    

            chr = [];
            mvx = [];
            wdur = [];
            man = [];
            hbd = [];
            hbm = [];
            baa = [];
            bap = [];
            haa = [];
            
            for s = 1:numel(ys),
                for w = wper{s}.data',
                    wdur(end+1) = diff(w);        
                    if strcmp(wper{1}.label,'walk'),
                        w = round(mean(w));
                    else,
                        w = w(1);
                    end

% $$$         chr(end+1,:) = (abs(ys{s}(w,:,1,2)))./(sqrt(ys{s}(w,:,1,1).*ys{s}(w,:,2,2)));
                    chr(end+1,:) = (abs(ys{s}(w,:,fetInds(1),fetInds(2))))./...
                        (sqrt(ys{s}(w,:,fetInds(1),fetInds(1)).*ys{s}(w,:,fetInds(2),fetInds(2))));
% $$$         baa(end+1,:) = angle(mean(ys{s}(w,:,:,1),3));
% $$$         bap(end+1,:) = mean(ys{s}(w,:,:,1),3).*conj(mean(ys{s}(w,:,:,1),3));        
% $$$         haa(end+1,:) = angle(mean(ys{s}(w,:,:,2),3));
                    hbd(end+1,:) = angle(ys{s}(w,:,fetInds(1),fetInds(2)));
                    hbm(end+1,:) = angle(mean(ys{s}(w,1:7,fetInds(1),fetInds(2)),2));
                    mvx(end+1,:) = (vxy{s}(w,:));
                    %man(end+1,:) = (ang{s}(w,5,7,2));
                    %man(end+1,:) = (fetp{s}(w,:));
                    %chr(end+1,:) = nanmean(abs(ys{s}(w,:,1,2)))./nanmean(sqrt(ys{s}(w,:,1,1).*ys{s}(w,:,2,2)));
                    %mvx(end+1,:) = nanmean(vxy{s}(w,:));
                end
            end

% $$$ figure,hist(chr(:,7),100)
% $$$  figure,hist(hbm(:),50)
% $$$ figure,imagesc(ts{1},fs{1},log10(ys{6}(:,:,1,1)')),axis xy; colormap jet;

            ind = wdur>5;
            out = [];
            mout = [];
            haout = [];
            baout = [];
            for f = 1:numel(fs{1}),
                out(:,f) = histc(chr(ind,f),linspace(0,1,60));
                mout(:,f) = histcounts(hbd(ind,f),30,'BinLimits',[-pi,pi]);    
% $$$     baout(:,f) = histcounts(baa(ind,f),30,'BinLimits',[-pi,pi]);        
% $$$     haout(:,f) = histcounts(haa(ind,f),30,'BinLimits',[-pi,pi]);    
            end

% $$$ figure();
% $$$ subplot(211);
% $$$ imagesc(linspace(-pi,pi,30),fs{1},haout')
% $$$ subplot(212);
% $$$ imagesc(linspace(-pi,pi,30),fs{1},baout')


% $$$ figure,imagesc(1:30,fs{1},out')

            hfig = figure();
            hfig.Units = 'centimeters';
            hfig.Position = [1,1,7,16];
            hfig.PaperPositionMode = 'auto';
            % PLOT phaseDiff vs frequency pdf
            hax = axes;
            hax.Units = 'centimeters';
            [~,hcb] = imagescnan({linspace(-pi,pi,30),fs{1},bsxfun(@rdivide,mout,sum(mout))'},[0,0.25],'linear',true);
            title({[fett{1}{fetInds(1)}],[ ' - ' fett{1}{fetInds(2)}]});
            hax.Position = [2,11,3,3];
            hcb.Units = 'centimeters';
            hcb.Position = [5.1,11,0.5,3];
            % PLOT histogram (solid bars) of mean phaseDiff distribution between 1-4Hz 
            hax = axes;
            hax.Units = 'centimeters';
            hax.Position = [2,7,3,3];
            histogram(hbm(ind),30,'BinLimits',[-pi,pi],'Normalization','probability');    
            title(wper{1}.label);
            ylim([0,0.25]);
            xlim([-pi,pi]);
            % PLOT histogram (stairs) of mean phaseDiff distribution between 1-4Hz 
            hax = axes;
            hax.Units = 'centimeters';
            hax.Position = [2,2,3,3];
            histogram(hbm(ind),30,'BinLimits',[-pi,pi],'DisplayStyle','stairs','Normalization','probability');
            title(hbfet{1}.label);
            ylim([0,0.25]);
            xlim([-pi,pi]);
            FigName = ['head_body-',hbfet{1}.label,'-',sprintf('%i_%i',fetInds),'-',wper{1}.label];
            print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
            print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));

            pause(1);
            delete(hfig);
        end
    end
end

figure,hold on
plot(sq(median(abs(cell2mat(cf(@(f,s) (f.segs(s{'n'}(:,1)-60,120,nan)), hbfet,stca))),2,'omitnan')))

plot(sq(median(abs(cell2mat(cf(@(f,s) (f.segs(s{'n'}(:,1)-60,120,nan)), hbfet,stc))),2,'omitnan')))


figure,
plot(sq(median(cell2mat(cf(@(f,s) (f.segs(s{'w'}(:,1)-60,120,nan)), hbfet,stca)),2,'omitnan')))


figure,plot(chr(ind,5),mvx(ind,2),'.');
figure,plot(chr(ind,5),man(ind,1),'.');

figure,plot(hbd(ind,5),man(ind,1),'.');


figure,hist2([chr(ind,3),mvx(ind,1)],20,20);
figure,hist2([chr(ind,3),man(ind,2)],10,10);

figure,hist2([chr(ind,3),man(ind,2)],10,10);





i = 4;
figure,
subplot(i,1,1);hold('on');imagesc(ts{1},fs{1},log10(ys{1}(:,:,3,3)')),axis('xy');colormap('jet');ylim([0,20]);plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,hbfet{1}(:,3)*100+10),
subplot(i,1,2);hold('on');imagesc(ts{1},fs{1},log10(ys{1}(:,:,2,2)')),axis('xy');colormap('jet');ylim([0,20]);plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,hbfet{1}(:,2)*100+10),
subplot(i,1,3);hold('on');imagesc(ts{1},fs{1},log10(ys{1}(:,:,1,1)')),axis('xy');colormap('jet');ylim([0,20]);plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,hbfet{1}(:,1)*100+10),plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,features{1}(:,26)+10)
subplot(i,1,4);hold('on');plotSTC(stc{1},1);ylim([0.5,7]);
linkaxes(findobj(gcf,'Type','axes'),'x');


%% STEP lateral sway analysis
OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_3/parts';



lspecParm = repmat({struct('nFFT',2^9,'Fs',features{1}.sampleRate,...
                           'WinLength',2^8,'nOverlap',2^8-8,...
                           'FreqRange',[1,6])},...
                   [1,numel(Trials)]);

[yl,fl,tl] = cf(@(t,f,p)  fet_spec(t,f,'mtcsdglong',true,'defspec',p,'overwrite',false),  Trials,lfet,lspecParm);


fxyz     = cf(@(x)    x.copy(),  xyz);
           cf(@(x)    x.filter('ButFilter',3,2.4,'low'),  fxyz);
vxySteps =            cf(@(x)    x.vel([],[1,2]),  fxyz);
           cf(@(x,y)  x.resample(y),           fxyz, yl);
vxy      = cf(@(x)    x.vel([],[1,2]),  fxyz);


i = 4;
figure,
subplot(i,1,1);hold('on');imagesc(ts{1},fl{1},log10(yl{1}(:,:,5,5)'));plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,lfet{1}(:,5)/5+4);
subplot(i,1,2);hold('on');imagesc(ts{1},fl{1},log10(yl{1}(:,:,4,4)'));plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,lfet{1}(:,4)/5+4);
subplot(i,1,3);hold('on');imagesc(ts{1},fl{1},log10(yl{1}(:,:,1,1)'));plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,lfet{1}(:,1)/5+4);
ForAllSubplots('caxis([-3,1.5]);ylim([1,8]);axis xy;colormap jet;');
subplot(i,1,4);hold('on');plotSTC(stc{1},1);ylim([0.5,7]);
linkaxes(findobj(gcf,'Type','axes'),'x');

nyl = cf(@(y) y.copy(), yl);
      cf(@(y) y.unity('drpOutPrctile',[2,98]), nyl);


i = 4;
figure,
subplot(i,1,1);hold('on');imagesc(ts{1},fl{1},(nyl{1}(:,:,5,5)'));plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,lfet{1}(:,5)/5+4);
subplot(i,1,2);hold('on');imagesc(ts{1},fl{1},(nyl{1}(:,:,4,4)'));plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,lfet{1}(:,4)/5+4);
subplot(i,1,3);hold('on');imagesc(ts{1},fl{1},(nyl{1}(:,:,1,1)'));plot((1:size(hbfet{1},1))./hbfet{1}.sampleRate,lfet{1}(:,1)/5+4);
ForAllSubplots('caxis([-2,2.5]);ylim([1,8]);axis xy;colormap jet;');
subplot(i,1,4);hold('on');plotSTC(stc{1},1);ylim([0.5,7]);
linkaxes(findobj(gcf,'Type','axes'),'x');
      
      

state = cf(@(s,y) s{'walk&gper'},   stc  ,nyl);
        cf(@(s,y) resample(s,y),    state,nyl);
%syl   = cf(@(y,s) log10(y(s.data(diff(s.data,1,2)>14,:),:,1,1)),   nyl ,state);
syl   = cf(@(y,s) (y(s.data(diff(s.data,1,2)>14,:),:,1,1)),   nyl ,state);
svel  = cf(@(v,s) v(s.data(diff(s.data,1,2)>14,:),'spine_lower'),  vxy,state);


state = cf(@(s,y) s{'walk&gper'},   stc  ,yl);
        cf(@(s,y) resample(s,y),    state,yl);
%syl   = cf(@(y,s) log10(y(s.data(diff(s.data,1,2)>14,:),:,1,1)),   nyl ,state);
syl   = cf(@(y,s) (y(s.data(diff(s.data,1,2)>14,:),:,1,1)),   yl ,state);
svel  = cf(@(v,s) v(s.data(diff(s.data,1,2)>14,:),'spine_lower'),  vxy,state);

syl   = cat(1,syl{:});
svel  = cat(1,svel{:});

%syl = unity(syl(nniz(syl),:));

% $$$ svel(svel<1e-3) = nan;
% $$$ svel = log10(svel);
svel(~nniz(svel))=nan;


vedgs = linspace(5,60,25);
%vedgs = linspace(0.5,1.75,50);
vbs = discretize(svel,vedgs,'IncludedEdge','right');
mrv = nan(numel(vedgs)-1,nyl{1}.size(2));
srv = nan(numel(vedgs)-1,nyl{1}.size(2));
% $$$ for f =1:size(syl,2),
% $$$     mrv(:,f) = accumarray(vbs(nniz(vbs)),syl(vbs(nniz(vbs)),f),[numel(vedgs),1],@nanmean);
% $$$ end

for f = 1:numel(vedgs)-1,
    mrv(f,:) = mean(syl(vbs==f,:),'omitnan');
    srv(f,:) = std(syl(vbs==f,:),'omitnan');
end



% PLOT lateral lower spine sway mean power spectrum as a function of body speed
figure();
hax = subplot2(4,1,1:3,1);
imagesc(vedgs,fl{1},(mrv)');
axis('xy');
colorbar();
%colormap('jet');
%caxis([0,0.6]);



% DETECT steps
flfet = cf(@(f) f.copy(),  lfet);
        cf(@(f) f.filter('ButFilter',3,[1.2,6],'bandpass'), flfet);
[steps,swayMagnitude] = cf(@(f) LocalMinima(-abs(f(:,1)),5,-3),  flfet);


for s = 1:numel(Trials),steps{s}(diff(steps{s})>90) = [];end


stepState = cf(@(s) s{'walk&gper'},   stc);
steps = cf(@(s,p) s(WithinRanges(s,p.data(diff(p.data,1,2)>14,:))),steps,stepState);

stepVels = cf(@(v,s) v(s,1) ,vxySteps,steps);

stepTimes = cat(1,steps{1});
stepSpeeds = cat(1,stepVels{1});

%ssl = log10(stepSpeeds);
ssl = (stepSpeeds);
ssf = 1./(([diff(stepTimes);1]+[1;diff(stepTimes)])).*120;

ssl = ssl(ssf>1);
ssf = ssf(ssf>1);

%figure,plot(ssl,ssf,'.')

vbins = linspace([vedgs([1,end]),10]);
for v = 1:numel(vbins)-1
    sslInd = ssl>vbins(v)&ssl<vbins(v+1);
    mssf(v) = mean(ssf(sslInd),'omitnan');
    sessf(v) = std(ssf(sslInd),'omitnan')./sqrt(sum(sslInd));
end

hold('on');
haxe = errorbar(vbins(1:end-1)+diff(vbins)/2,mssf,sessf);
haxe.Color = [0,0,0];
haxe.LineWidth = 1;

%imhand = imagesc(vedgs,fs,bsxfun(@plus,mrv,abs(min(mrv(:))))');


% PLOT marginal distribution of log10 speed
haxd = subplot2(4,1,4,1);
hist(svel,100);
haxd.Position(3) = hax.Position(3);
xlim([vedgs([1,end])])

OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'Shared/Behavior Paper/Figures/Figure_3/parts';

FigName = ['lateral_sway_freq_pow'];
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(OwnDir,FigDir,[FigName,'.png']));
