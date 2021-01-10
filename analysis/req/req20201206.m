;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.

figure();
plot(xyzp(tper,1),xyzp(tper,2),'.');
Lines([],[-500:100:500],'k');
Lines([-500:100:500],[],'k');

%req20201117(Trial)

% Interneuron summary
figure()
for u = unitsInt,
    
end


MjgER2016_load_data();

phzCorrection = [pi*1.25,pi*1.25,                             ... er01
                 pi/2,pi/2,pi/2,                    ... ER06
                 pi/1.25,pi/1.25,                             ... Ed10
                 0,0,0,0,0,0,0,0,0,                 ... jg04                 
                 pi/4,pi/4,pi/4,pi/4,pi/4,pi/4,pi/4 ... jg05
                 pi/4,pi/4,pi,pi/1.25,pi*1.25]; % new units - jg05, jg05, ER06, Ed10, er01

cf(@(t,p) report_interneuron_summary(t,'phzCorrection',p), Trials(end), num2cell(phzCorrection(end)));


unitsPyr = {...
    [161],...           er01 20110719
    [],...              er01 20110721
    [],...              ER06 20130612
    [],...              ER06 20130613
    [115,118,152],...   ER06 20130614
    [],...              Ed10 20140816
    [],...              Ed10 20140817
    [],...              jg04 20120128
    [1],...             jg04 20120129
    [18],...            jg04 20120130
    [22],...            jg04 20120131
    

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
    

t = 20;    
Trial = Trials{t};
spk = Trial.spk.copy();
spk.create(Trial,              ...  
           sampleRate,         ...
           '',                 ...
           unitsInts{t},       ...
           ''                  ...
);    
pft = pfs_2d_theta(Trial);

stc = Trial.stc.copy();    
rper = [stc{'R&s',sampleRate}];
sper = [stc{'s',sampleRate}];
rper = [stc{'r',sampleRate}];
rper.data(diff(rper.data,1,2)<200,:) = [];
%sper = [Trial.stc{'t',sampleRate}];
rper = stc.get_state_transitions(Trial,{'rear','pause'},0.2,xyz);
rper = rper-80;
%pch = fet_HB_pitchB(Trial,sampleRate);

figure();
plot(xyz(:,'spine_upper',3));
Lines(rper(:),[],'k');

frq = copy(phz);
frq.data = circ_dist(frq.data,circshift(frq.data,1))/(2*pi)*sampleRate;    
figure,plot(frq.data)


    


figure,
for u = unitsInts{t};
res = spk(u);    
res = res(WithinRanges(res,rper));
clf();
scatter(frq(res),phz(res),10,xyz(res,'spine_upper',3),'filled');
colormap('jet');    
xlim([0,15]);
waitforbuttonpress();
end

figure();
for u = unitsInts{t};
res = spk(u);
%res = res(WithinRanges(res,sper));
%[mccg,tbins] = CCG([mean(rper.data,2);res],[ones([size(rper,1),1]);2*ones(size(res))],1,100,sampleRate,[1,2]);
[mccg,tbins] = CCG([rper(:,1);res],[ones([size(rper,1),1]);2*ones(size(res))],4,70,sampleRate,[1,2]);
subplot(2,2,[1]);plot(pft,u);
subplot(2,2,[2]);bar(tbins,mccg(:,1,2));
[mccg,tbins] = CCG([rper(:,2);res],[ones([size(rper,1),1]);2*ones(size(res))],4,70,sampleRate,[1,2]);
subplot(2,2,[4]);bar(tbins,mccg(:,1,2));
title(num2str(u));
waitforbuttonpress();
end


sper = [Trial.stc{'t-m-s',sampleRate}];
sper = [Trial.stc{'s-t',sampleRate}];
sper = [Trial.stc{'s-t',sampleRate}];

sper = [Trial.stc{'R-t',sampleRate}]+[-0.2,0.2];
u = [ 92,106];
figure();
for x = 1:14,
    for y = 1:14,
        u = [ unitsInts{t}(x),unitsInts{t}(y)];
        res1 = spk(u(1));res1 = res1(WithinRanges(res1,sper));
        res2 = spk(u(2));res2 = res2(WithinRanges(res2,sper));
        [mccg,tbins] = CCG([res1;res2],[ones(size(res1));2*ones(size(res2))],1,70,sampleRate,[1,2]);
        subplot2(14,14,y,x);
        bar(tbins,mccg(:,1,2));
        if x == y,
            ylim([0,75]);
        end
    end
end



figure,
for u = unitsInts{t};
res = spk(u);    
res = res(WithinRanges(res,sper));
pshift = -8*pi:pi/4:8*pi;
%pshift = -250:25:250;
mccg = zeros([101,numel(pshift)]);
for p = 1:numel(pshift)
rper = LocalMinima(circ_dist(phz.data,pshift(p))+pi,10,0.1);
%rper = LocalMinima(circshift(phz.data,pshift(p))+pi,10,0.1);
[tccg,tbins] = CCG([rper(:,1);res],[ones([size(rper,1),1]);2*ones(size(res))],4,50,sampleRate,[1,2]);
mccg(:,p) = tccg(:,1,2);
end;
clf();
imagesc(mccg');
Lines(51,[],'k');
Lines([],numel(pshift)/2,'k');
colormap('jet');
axis('xy');
waitforbuttonpress();
end


figure,
for u = unitsInts{t};
res = spk(u);    
res = res(WithinRanges(res,sper));
tic
rper = LocalMinima(circ_dist(phz.data,-circ_mean(phz(res))),10,0.1);
toc
[mccg,tbins] = CCG([rper(:,1);res],[ones([size(rper,1),1]);2*ones(size(res))],4,50,sampleRate,[1,2]);
clf();
subplot(121);
rose(phz(res),36);
subplot(122);hold('on');
bar(tbins,mccg(:,1,2)-mean(mccg(:,1,2)));
plot(tbins,normalize(RectFilter(conv((mccg(:,1,2)-mean(mccg(:,1,2))).^2,ones([11,1]),'same'),3,5),'range'),'r');
Lines(tbins(51),[],'k');
waitforbuttonpress();
end



    
cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials(1),unitsInts(1));
tag = 'interneurons_xyhb_2020';
tag = 'interneurons_xyhb_2020_bhv'; overwrite = true;
pfi = cf(@(t,u) req20201117(t,u,tag,false,false,overwrite), Trials,unitsInts);



[rmaps,cluSessionMap] = decapsulate_and_concatenate_mtaapfs(pfi,unitsInts);

MjgER2016_load_data();
configure_default_args();

% LOAD place restricted behavior fields
bfs   = cf(@(t,u)   compute_bhv_ratemaps(t,u),  Trials, units);
%pfsa = cf(@(s,t) cat(2,{t},s), pfss,pfts);
% COMPUTE bfs erpPCA
[eigVecs, eigScrs, eigVars, unitSubset, validDims, zrmMean, zrmStd] = ...
                    compute_bhv_ratemaps_erpPCA(bfs, units);
numComp = size(eigVecs,2);
fpc  = cell([1,numComp]);
for i = 1:numComp,
    fpc{i} = nan(size(validDims));
    fpc{i}(validDims) = eigVecs(:,i);
end
fpcLims = [min(cellfun(@min,fpc)),max(cellfun(@max,fpc))];

figure,
imagesc(reshape_eigen_vector(fpc{4},bfs(1)))
axis('xy');



bhvMask = false(size(validDims));
bhvMask(validDims) = true;
bhvMask = reshape_eigen_vector(bhvMask,bfs)';
bhvLims = [-1.6, 0.6; ...
           -0.5, 1.7];




whos rmaps

bins = pfi{1}.adata.bins;
binDims = pfi{1}.adata.binSizes';

rmapa = nan([binDims,size(rmaps,2)]);
for u = 1:size(rmaps,2),
    rmapa(:,:,:,:,u) = reshape(rmaps(:,u),binDims);
end


rmap = rmapa(:,:,:,:,200);
figure();
for x = 1:6,
    for y = 1:6,
        subplot2(6,6,7-y,x);
        imagescnan({bins{3:4},sq(rmap(x,y,:,:))'});
        axis('xy');
    end
end

bmaps = reshape(rmapa(3,3,:,:,:),[],size(rmapa,length(binDims)+1));

figure,
bar(sum(double(~isnan(bmaps)),2));
% Valid Dimensions
disp(sum(double(sum(double(~isnan(bmaps)),2)>100)))

validDims = sum(double(~isnan(bmaps)),2)>200;

figure();
imagesc(reshape(validDims,binDims(end-1:end))');
axis('xy');


vmaps = bmaps(validDims,:);
vmaps(:,sum(isnan(vmaps))>5)= [];
vmaps(isnan(vmaps(:)))=0;

[LU,LR,FSr,VT] = erpPCA(vmaps',5);

figure()
for v = 1:5,
    eigVec = nan(binDims(end-1:end));
    eigVec(validDims) = LU(:,v);
    subplot(1,5,v);
    imagescnan({bins{end-1:end},eigVec'});
    axis('xy');
end



sper = [Trial.stc{'t-m-s',sampleRate}];

thpks = LocalMinima(abs(phz.data),10,0.1);
thpks = thpks(WithinRanges(thpks(:,1),sper),:);
thpks = [thpks,circshift(thpks,-1)];
thpks([1,end],:) = [];
thpks(diff(thpks,1,2)>50,:) = [];

%figure(); hist(diff(thpks,1,2),1000)
u = 76;
res = spk(u);
thspkPPC = nan([size(thpks,1),1]);
thspkN = nan([size(thpks,1),1]);
thspkCM = nan([size(thpks,1),1]);
for t = 1:size(thpks,1),
    tres = res(WithinRanges(res,thpks(t,:)));
    thspkN(t) = numel(tres);
    if thspkN(t)>0,
        thspkCM(t)  = circ_mean(phz(tres));
        thspkPPC(t) = PPC(phz(tres));
    end
end


figure();
ind = ~isnan(thspkPPC)&~isnan(thspkN);
hist2([thspkPPC(ind),thspkN(ind)],30,15);

figure();
ind = ~isnan(thspkPPC) & ~isnan(thspkN) & thspkN>3;
out = hist2([thspkPPC(ind),thspkCM(ind)],30,15);
imagesc([out,out]')

figure();
ind = ~isnan(thspkPPC) & ~isnan(thspkN) & thspkN>2;
out = hist2([thspkN(ind),thspkCM(ind)],15,30);
imagesc([out,out]')

dthpks = diff(thpks,1,2);

figure();
for s  = -2:2,
    subplot(1,5,3+s);
    tdthpks = circshift(dthpks,s);
    ind = ~isnan(thspkPPC) & ~isnan(thspkN) & thspkN>=2;
    out = hist2([tdthpks(ind),thspkCM(ind)],40,60);
    imagesc([out,out]');
    caxis([0,40]);
    title(num2str(u));
end
colormap('jet');


%hist2([diff(thpks,1,2),thspkPPC'],30,30);



% fet

function [fet,featureTitles,featureDesc] = fet_HB_pitchB(Trial,varargin)
% function [fet,featureTitles,featureDesc] = fet_head_pitch(Trial,varargin)
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%     referenceTrial'  , 'Ed05-20140529.ont.all'
%     referenceFeature', ''
%

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate'   , Trial.xyz.sampleRate,                                       ...
                 'normalize'       , false,                                                      ...
                 'procOpts'        , {{}},                                                       ...
                 'referenceTrial'  , '',                                                         ...
                 'referenceFeature', '',                                                         ...
                 'overwriteFlag'   , false                                                       ...
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature,overwriteFlag] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'Head body pitch','fet_HB_pitch','b');                  

if overwriteFlag,
% XYZ preprocessed 
    xyz = Trial.load('xyz');
% LOAD pitches 
    pch = fet_HB_pitch(Trial);
    if ~isempty(referenceTrial),
% MAP to reference trial
        pch.map_to_reference_session(Trial,referenceTrial,referenceFeature);
    end
    pch.resample(newSampleRate);
    xyz.resample(newSampleRate);
% CONCATENATE features
    fet.data = [circ_dist(pch(:,3),pch(:,1)),pch(:,1)];
    fet.data(~nniz(xyz),:)=0;
    fet.save();
else
    load(fet,Trial);
end


featureTitles = {};
featureDesc = {};
if nargout>1,
    featureTitles(end+1) = {'Pitch HCHN-BMBU'};    
    featureDesc(end+1) = {['head pitch relative to xy plane']};
    featureTitles(end+1) = {'Pitch BMBU'};    
    featureDesc(end+1) = {['upper body pitch relative to xy plane']};
end

% END MAIN -----------------------------------------------------------------------------------------