
MTAConfiguration('/gpfs01/sirota/bach/data/gravio/','absolute');

%% place field fig

Trial = MTATrial('jg05-20120317');
%Trial = MTATrial('er06-20130614','all-cof');
if Trial.xyz.isempty,Trial.xyz.load(Trial);end


%units = [3:30];
units = 1:95;
states = {'rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
numsts = numel(states)

xyo={};pfs={};accg={};
for i = 1:numsts,
    %For MTAAknnpf
    pfs{i}  =        MTAAknnpfs(Trial,units,states{i},0,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
    %[xyo{i},xyb]    = xytrajocc(Trial,Trial.xyz(Trial.stc{states{i}},Trial.trackingMarker,[1,2]),pfs{i}.parameters.binDims,pfs{i}.parameters.distThreshold,'xy');
    %[accg{i},tbin] = autoccg   (Trial,units,states{i});


    %For MTAApfs
    %pfs{i}  =        MTAApfs(Trial,units,states{i},0,'numIter',100);
    %[xyo{i},xyb]    = xytrajocc(Trial,Trial.xyz(Trial.stc{states{i}},Trial.trackingMarker,[1,2]),[20,20],70,'xy');
    %[accg{i},tbin] = autoccg   (Trial,units,states{i});
end

for u = units
for state = 1:numsts
    pfstats(state,u) = PlaceFieldStats(Trial,pfs{state},u);    
end
end

for u = units
    %figure,
    %u = 15;

for state = 1:numsts

    % xy occupancy
    subplot2(6,numsts,1,state);
    imagescnan({xyb{1},xyb{2},xyo{state}},[],0,1);

    subplot2(6,numsts,2,state);
    bar(tbin,accg{state}(:,u));axis tight
    
    subplot2(6,numsts,3,state);
    pfs{state}.plot(u,[],[],[0,max([pfstats(:,u).peakFR])])
    
    subplot2(6,numsts,4,state);
    pfs{state}.plot(u,'sig');

    subplot2(6,numsts,5,state);
    hist(pfstats(state,u).shuffledPatchArea,100);
    Lines(pfstats(state,u).shuffledPatchArea(1),[],'r');
end

saveas(gcf,fullfile('/gpfs01/sirota/bach/homes/gravio/figures/knnpf',[Trial.filebase '.knnpf-' num2str(u) '.png']),'png')

%set(gcf,'Position',get(2,'Position'));
%set(gcf,'Name',['unit ' num2str(u)])
end


%| xyocc_a  |  xyocc_t  |  xyocc_r&t  | xyocc_w&t    | xyocc_h&t    |  xyocc_l&t    |
%| accg_all |  accg_t   |  accg_r&t   | accg_w&t     | accg_h&t     |  accg_l&t     |
%| pf_all   |  pf_theta |  pf_rear&th | pf_walk&th   | pf_hwalk&th  |  pf_lwalk&th  |
%| pf_a_sig |  pf_t_sig |  pf_r&t_sig | pf_w&t_sig   | pf_h&t_sig   |  pf_l&t_sig   |
%| phst_a   |  phst_t   |  phst_r&t   | phst_w&t     | phst_h&t     |  phst_l&t     |


states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};

for state = 1:numsts
Trial.spk.create(Trial,Trial.sampleRate,states{state});
%[tccg,t] = Trains2CCG({Trial.spk.res},{Trial.spk.clu},16,60,Trial.sampleRate,'count');
[tccg,t,pairs] = CCG(Trial.spk.res,Trial.spk.clu,16,60,Trial.sampleRate,Trial.spk.map(:,1),'count',[]);
figure

us1 = 18:25;
us2 = 18:25;
icount = 1;
for i = us1
jcount = 1;
for j = us2
subplot2(numel(us1),numel(us2),icount,jcount);
bar(t,tccg(:,j,i)),axis tight,
title([num2str(j) ' _ ' num2str(i)]);

jcount = jcount + 1;
end
icount = icount + 1;
end
set(gcf,'Position',get(15,'Position'));
set(gcf,'Name',states{state})
end



reshape(cat(1,pfstats(:,:).patchCOM),[5,numel(units),2])

figure,plot(reshape[pfstats([2,3],:).peakFR],2,[])','.')


figure,plot(log10([pfstats([3],:).peakFR]),log10([pfstats([2],:).peakFR]),'.')
line([0;2],[0;2])

figure,plot(log10([pfstats([4],:).peakFR]),log10([pfstats([5],:).peakFR]),'.')
line([0;2],[0;2])



%% peak firing rate diff vs patch distance
u = 20;
%sid = ~cellfun(@isempty,regexp(states,'(hwalk&theta)|(lwalk&theta)'));
sid = ~cellfun(@isempty,regexp(states,'(^rear&theta?)|(^walk&theta?)'));
figure
hold on
for u = units([Trial.nq.SpkWidthR(1:29)]>.5&[Trial.nq.eDist(units)]>30),
plot(pdist(cat(1,pfstats(sid,u).patchCOM)),diff([pfstats(sid,u).patchPFR]),'.')
end



%% RateMap Correlations
Trial = MTATrial('jg05-20120317');
units = 1:size(Trial.spk.map,1);
states = {'rear&theta','walk&theta','hwalk&theta','lwalk&theta'};
numsts = numel(states);

pfs={};
for i = 1:numsts,
    pfs{i}  =        MTAAknnpfs(Trial,units,states{i},0,'numIter',1000,'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
end


%prinunits = Trial.selectUnits({{Trial.nq,'SpkWidth
stsCor = zeros(numel(units),numsts,numsts,2,2,2);

for u = units,
    for statei = 1:numsts,
        for statej = 1:numsts,
            
            hsig = pfs{statei}.plot(u,'sig');
            lsig = pfs{statej}.plot(u,'sig');
            hmap = pfs{statei}.plot(u);
            lmap = pfs{statej}.plot(u);
            
            mind = hsig<.05&lsig<.05;
            if sum(mind(:))>10,
                [stsCor(u,statei,statej,:,:,1),stsCor(u,statei,statej,:,:,2)] = corrcoef(hmap(mind),lmap(mind));
            else
                stsCor(u,statei,statej,:,:,1) = zeros(2,2);
                stsCor(u,statei,statej,:,:,2) = ones(2,2);
            end
        end
    end
end



cind = stsCor(:,2,3,1,2,2)<0.05&stsCor(:,2,4,1,2,2)<0.05;
figure,plot(stsCor(cind,2,3,1,2,1),stsCor(cind,2,4,1,2,1),'.')
xlim([-1,1]),ylim([-1,1])
line([-1;1],[-1;1])


cind =stsCor(:,1,2,1,2,2)<0.05&stsCor(:,1,3,1,2,2)<0.05;
figure,plot(stsCor(cind,1,2,1,2,1),stsCor(cind,1,3,1,2,1),'.')

figure,plot(sq(max(pfs{1}.data.rateMap(:,:,1))),sq(max(pfs{2}.data.rateMap(:,:,1))),'.')
figure,plot(log10(sq(max(pfs{3}.data.rateMap(:,:,1)))),log10(sq(max(pfs{4}.data.rateMap(:,:,1)))),'.')

figure,plot(log10(sq(max(pfs{3}.data.rateMap(:,:,1)))),log10(sq(max(pfs{4}.data.rateMap(:,:,1)))),'.')


u=29;
s1 = 1;
s2 = 2;

for u = units,
hsig = pfs{s1}.plot(u,'sig');
lsig = pfs{s2}.plot(u,'sig');
hmap = pfs{s1}.plot(u);
lmap = pfs{s2}.plot(u);
mind = hsig<.05&lsig<.05;
%figure,
subplot(131);
plot(hmap(mind)/max(hmap(:)),lmap(mind)./max(lmap(:)),'.');
xlabel(pfs{s1}.parameters.states);
ylabel(pfs{s2}.parameters.states);
title('Rate Map Correlation')
xlim([0,1]),ylim([0,1.0])
line([0;1],[0;1])
subplot(132);
pfs{1}.plot(u);
title(pfs{s1}.parameters.states);
subplot(133);
pfs{2}.plot(u);
title(pfs{s2}.parameters.states);
saveas(gcf,fullfile('C:\Users\justi_000\Dropbox\figures\rmcorr',[Trial.filebase '.rmcorr_' pfs{s1}.parameters.states(~ismember(pfs{s1}.parameters.states,'&')) 'X' pfs{s2}.parameters.states(~ismember(pfs{s2}.parameters.states,'&')) '-' num2str(u) '.png']),'png')
end

