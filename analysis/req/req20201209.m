Trial = MTATrial.validate('jg05-20120312.cof.all');
stc = Trial.load('stc','msnn_ppsvd_raux');
phzCorrection = pi/4;
thetaRef = 69;

unitsInts = {...
    [ 31, 78, 82,125,169,195,203],...                                    er01 20110719
    [ 31, 76,105,147],...                                                er01 20110721
    [ 27, 32, 68, 69,105,124,125,126,128,182,220,221,222,225],...        ER06 20130612
    [ 10, 13, 15, 18, 60, 93,112,113,192,213,215,216,220],...            ER06 20130613
    [  5, 10, 11, 12, 22, 43, 49, 91, 94, 95,110,124,144,145,146,155,... ER06 20130614
       156,168,169,170,171,172,176,177,179,181,182],...                  
    [  9, 11, 12, 29, 30, 43, 45, 46, 51, 52, 54, 55, 56, 87],...        Ed10 20140816
    [ 11, 12, 31, 47, 48, 51, 56, 80, 81, 89],...                        Ed10 20140817
    [  8,  9, 16],...                                                    jg04 20120128
    [ 21],...                                                            jg04 20120129
    [ 24],...                                                            jg04 20120130
    [ 10, 24, 27],...                                                    jg04 20120131
    [ 10],...                                                            jg04 20120201
    [  4,  5],...                                                        jg04 20120210
    [  4,  5],...                                                        jg04 20120211
    [  6],...                                                            jg04 20120212
    [  2],...                                                            jg04 20120213
    [  5, 10, 15, 27, 28, 38, 64, 66,100,114,116,117,121,122],...        jg05 20120309
    [  4,  6,  7,  8, 27, 28, 43, 59, 71, 86, 99,100,101,102],...        jg05 20120310    
    [  2,  7,  9, 25, 41, 44, 46, 70, 71, 72, 98,112,113,118,...         jg05 20120311
     119,135,144,148,188,193,203,204,205],...
    [  3,  7,  8, 15, 16, 43, 45, 50, 76, 77, 92,106,124,184],...        jg05 20120312
    [  1,  5, 34, 60, 69],...                                            jg05 20120315
    [ 17, 28, 49, 55],...                                                jg05 20120316
    [ 15, 16, 17, 52],...                                                jg05 20120317
    [  2,  4, 20, 30, 39, 40, 51],...                                    jg05 20120323
    [  3,  7, 14, 16, 34],...                                            jg05 20120324
    [  5, 13, 27, 66, 71, 72, 75, 91, 94,105,127,157,209,232,233,251,254],...ER06 20130624
    [  4, 17, 20, 22, 23, 25, 29 ,30, 31, 40, 66, 69],...                Ed10 20140815
    [ 93,175] ...                                                         er01 20110722
};
    


% LOAD local field potential (lfp)
% RESAMPLE lfp to xyz sample rate
% COMPUTE lfp phase in theta band (6-12 Hz)
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',65:96);
phz = lfp.phase([5,13]);    

spk = Trial.spk.copy();
spk.create(Trial,              ...  
           lfp.sampleRate,         ...
           '',                 ...
           [],       ...
           ''                  ...
);    


sper = [Trial.stc{'t-m-s',lfp.sampleRate}];


%            1    2     3     4     5     6     7     8      9      10
frqBins = [1,4;5,12;14,22;26,32;32,46;46,60;60,74;74,88;88,112;112,126];
states = {'theta-groom-sit','sit&theta','sit-theta','rear&theta','hloc&theta','hpause&theta','lloc&theta',...
          'lpause&theta','pause-theta'};
cr = nan([size(lfp,2),size(frqBins,1),numel(states),size(spk.map,1)]);
cm = nan([size(lfp,2),size(frqBins,1),numel(states),size(spk.map,1)]);
for f = 1:size(frqBins,1),
    disp(frqBins(f,:))
    phz = lfp.phase(frqBins(f,:));    
    for s = 1:numel(states),
        sper = [Trial.stc{states{s},lfp.sampleRate}];
        for k = 1:size(spk.map,1),
            res = spk(k);
            res = res(WithinRanges(res,sper.data));
            if numel(res)>50,
                for c = 1:32,
                    cr(c,f,s,k) = circ_r(phz(res,c));
                    cm(c,f,s,k) = circ_mean(phz(res,c));
                end
            end
        end        
    end
end

numStates = numel(states);

figure();
for s = 1:numStates,
    subplot2(numStates,2,s,1);hold('on');
    for k = unitsInts{20},
        plot(cr(:,2,s,k),fliplr(1:32));
    end
    subplot2(numStates,2,s,2);hold('on');
    for k = unitsInts{20},
        plot(cm(:,2,s,k),fliplr(1:32));
    end
    xlim([-pi,pi]);
end

figure,
plot(sq(diff(cr([10,18],2,1,unitsInts{20}))),sq(cm(6,2,1,unitsInts{20})),'.');
hold('on');
plot(sq(diff(cr([10,18],2,1,92))),sq(cm(6,2,1,92)),'.');
plot(sq(diff(cr([10,18],2,1,45))),sq(cm(6,2,1,45)),'.');



figure();
for k = unitsInts{20};
    clf();
    xlr = 0.08;%max(max(max(sq(cr(:,:,:,k)))));
    if xlr>0 
        for s = 1:numStates,
            subplot2(numStates,2,s,1);hold('on');
                imagesc(bsxfun(@rdivide,cr(:,:,s,k),sum(cr(:,:,s,k),1)));
                caxis([0,xlr]);
                colormap(gca(),'jet');
                colorbar(gca());
                axis('ij');
                ylim([0,33]);
                title(states{s});
                Lines([],6 ,'k');                
                Lines([],10,'k');
                Lines([],16,'k');
                Lines([],21,'k');                
                Lines([],27,'k');     
                Lines(4.5,[],'k');
                Lines(6.5,[],'k');                
                Lines(7.5,[],'k');                
            subplot2(numStates,2,s,2);hold('on');
                imagesc(cm(:,:,s,k));
                caxis([-pi,pi]);
                colormap(gca(),'hsv')
                colorbar(gca());
                axis('ij'); 
                ylim([0,33]);
                Lines([],6 ,'k');
                Lines([],10,'k');
                Lines([],16,'k');
                Lines([],21,'k');                
                Lines([],27,'k');
                Lines(4.5,[],'k');
                Lines(6.5,[],'k');                
                Lines(7.5,[],'k');                
            subplot2(numStates,2,s,2);hold('on');
        end
    end
    title(num2str(k));
    waitforbuttonpress();
end


figure();
plot(sq(cm(:,2,1,unitsInts{20})),circ_dist(sq(cm(:,2,1,unitsInts{20})),sq(cm(:,2,3,unitsInts{20}))),'.');

figure();
plot(sq(cm(6,2,1,unitsInts{20})),mean(circ_dist(sq(cm(:,2,1,unitsInts{20})),sq(cm(:,2,3,unitsInts{20})))),'.');

figure();
plot(sq(cm(20,2,1,unitsInts{20})),(circ_dist(sq(cm(20,2,1,unitsInts{20})),sq(cm(20,2,3,unitsInts{20})))),'.');


figure();
plot(sq(cm(6,2,1,[45,92])),(circ_dist(sq(cm(20,2,1,[45,92])),sq(cm(20,2,3,[45,92])))),'.');


figure();
for k = 1:180;
    clf();
    xlr = max(max(sq(cr(:,f,:,k))));
    if xlr>0 
        for s = 1:numStates,
            subplot2(numStates,2,s,1);hold('on');
                plot(cr(:,f,s,k),fliplr(1:32)');
                xlim([0,xlr]);
            subplot2(numStates,2,s,2);hold('on');
                plot(cm(:,f,s,k),fliplr(1:32)');
            xlim([-pi,pi]);
        
        end
    end
title(num2str(k));
waitforbuttonpress();
end

% $$$ phz.data = unwrap(phz.data);
% $$$ phz.resample(xyz);    
% $$$ %phz.data = mod(phz.data+pi,2*pi)-pi;
% $$$ phz.data = mod(phz.data+pi,2*pi)-pi + phzCorrection; 
% $$$ phz.data(phz.data>pi) = phz.data(phz.data>pi)-2*pi;
% $$$ lfp.resample(xyz);


res = 
mtptchd(lfp.data,res,clu,2^8,lfp.sampleRate,2^7,2^6,

rlfp = Trial.load('lfp',61);
[ys,fs,ts] = fet_spec(Trial,rlfp,[],false,'defspec',struct('nFFT',2^7,'Fs',rlfp.sampleRate,...
                                                  'WinLength',2^6,'nOverlap',2^6*.875,...
                                                  'FreqRange',[100,300]));

figure();
imagesc(ts,fs,log10(ys.data)');
axis('xy');
tper = stc{'t'};
tper.cast('TimeSeries');
tper.resample
tys = ys.copy();
ysx = ys.copy();

xyz = preproc_xyz(Trial);
ysx.resample(xyz);
vxy = vel(xyz,{'spine_lower','hcom'},[1,2]);
ang = create(MTADang,Trial,xyz);

marker = 'spine_upper';

ind = log10(ysx(:,7))>2&xyz(:,marker,3)>0&vxy(:,1)>0.01;

figure,
%hist2([real(log10(ysx(ind,7))),xyz(ind,marker,3)],100,500)
%hist2([real(log10(ysx(ind,7))),log10(vxy(ind,1))],100,100)
fb = 1:3:20;
for f = 1:7
subplot(1,7,f)
%hist2([real(log10(ysx(ind,fb(f)))),ang(ind,5,7,2)],linspace(2,6,50),linspace(-pi/2,pi/2,100));
hist2([real(log10(ysx(ind,fb(f)))),log10(vxy(ind,2))],linspace(2,6,50),linspace(-2,2,100),'xprob')
caxis([0,0.008]);
%caxis([0,500]);
end


figure();
plot(log10(ysx(ind,7)),xyz(ind,'hcom',3),'.')





[mccg,tbins] = CCG(spk.res,spk.clu,1,200,spk.sampleRate,spk.map(:,1));


figure,
for u = 92:120,
    bar(tbins,mccg(:,u,92));
    title(num2str(u));
    waitforbuttonpress();
end

figure();
bar(tbins,mccg(:,106,85));
xlim([-50,50]);

pft = pfs_2d_theta(Trial);
pfs = pfs_2d_states(Trial,spk.map(:,1));

slabels ={'loc','lloc','hloc','rear','pause','lpause','hpause'};

figure();
for u = spk.map(37:end,1)';
    clf();
subplot(1,9,1);
bar(tbins,mccg(:,u,u));
title(num2str(spk.map(u,:),'u:%i el:%i clu:%i'));
axis('tight');
xlim([-50,50]);
subplot(1,9,2);
plot(pft,u,2,'text');
title(['theta']);
for s = 1:7,
    subplot(1,9,s+2);
    plot(pfs{s},u,1,'text');
    title([slabels{s}]);
end
waitforbuttonpress();
end
Trial.load('nq');

ur = [5,93,102,109,128];
[ur',pft.data.si(ur)',mrt(ur),Trial.nq.eDist(ur),Trial.nq.Refrac(ur),nq_type(ur),nq_quality(ur),nq_type(ur)<0,nq_quality(ur)>0,nq.Refrac(ur)&0.001]    


figure,plot(pft.data.si,mrt','.','MarkerSize',10)
hold('on');plot(pft.data.si(units{20}),mrt(units{20})','.r','MarkerSize',10)
hold('on');plot(pft.data.si(unitsInts{20}),mrt(unitsInts{20})','.g','MarkerSize',10)
hold('on');plot(pft.data.si(ur),mrt(ur)','.c','MarkerSize',10)

figure,
plot(pft.data.si,mrt','.','MarkerSize',10);
ind = units{20};    hold('on');plot(pft.data.si(ind),mrt(ind)','.r','MarkerSize',10);
ind = unitsInts{20};hold('on');plot(pft.data.si(ind),mrt(ind)','.g','MarkerSize',10);
ind = ur;           hold('on');plot(pft.data.si(ind),mrt(ind)','.c','MarkerSize',10)

figure(); hold('on');
ind = mrt>0.4;        plot(nq.SNR(ind),mrt(ind)','.','MarkerSize',10); Lines([],0.4,'k');
ind = units{20};      plot(nq.SNR(ind),mrt(ind)','.r','MarkerSize',10);
ind = unitsInts{20};  plot(nq.SNR(ind),mrt(ind)','.g','MarkerSize',10);
ind = ur;             plot(nq.SNR(ind),mrt(ind)','.c','MarkerSize',10)
ind = uc;             plot(nq.SNR(ind),mrt(ind)','.y','MarkerSize',10)



% select_placefields
uc = [38,39,54,55,56,82,94,112,121,125,143,146];
[uc',Trial.nq.eDist(uc),Trial.nq.Refrac(uc)]


Trial.load('nq');    
mrt = pft.maxRate(spk.map(:,1));

load(fullfile(MTASession([]).path.cfg,'unit_selection_criteria.mat'));    
nq = Trial.nq;
nq_type    = nq.(usp.type.fields{2})   -polyval(usp.type.pram,   nq.(usp.type.fields{1}));
nq_quality = nq.(usp.quality.fields{2})-polyval(usp.quality.pram,nq.(usp.quality.fields{1}));


    
[uc',Trial.nq.eDist(uc),Trial.nq.Refrac(uc),nq_type(uc),nq_quality(uc),nq_type(uc)<0,nq_quality(uc)>0,nq.Refrac(uc)&0.001]    

units = select_units(Trial,'pyr');
units = units();
spk.create(Trial,[],'theta-groom-sit',units,'deburst');
spkCnt = accumarray(spk.clu,ones([numel(spk.clu),1]),[size(spk.map,1),1],@sum);
units = units(spkCnt(units)>minSpkCnt);
%save(filename,'units');
