
Trial = MTATrial('jg05-20120317','all');
Trial.load('nq');



stc_mode = 'auto_wbhr';
sesList = {{'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};

states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};


%pftype = 'MTAAknnpf';

numsts = numel(states);
pfstats = {};
pfshuff = {};
%pfs={};

try,matlabpool('open',12);,catch,matlabpool('close');matlabpool('open',12),end

map = [];
anq = {};
for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
    Trial.stc.updateMode(stc_mode);Trial.stc.load;
    Trial.load('nq');
    %Trial.load('xyz');


    units = select_units(Trial,18,'pyr');
    map = cat(1,map,cat(2,size(map,1)+[1:numel(units)]',ses*ones([numel(units),1]),Trial.spk.map(units,:)));

    tnq = StructArray(Trial.nq,1);
    anq{ses} = tnq(units);

    clear tpfstats,
    for i = 1:numsts,
        pfs =     MTAAknnpfs(Trial,units,states{i},0,'numIter',1000, ...
                            'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
        clear ts,clear tss,
        parfor u = 1:numel(units),
            [ts(u)] = PlaceFieldStats(Trial,pfs,units(u));
        end      
        tpfstats{i} = CatStruct(ts,[],1,[],0);
    end

    pfstats{ses} = CatStruct(cat(1,tpfstats{:}),[],2,[],0);
end

matlabpool('close');

apfs = CatStruct(cat(1,pfstats{:}),[],1,[],0);
anq = CatStruct(cat(1,anq{:}),[],1,[],0)


% $$$ u  = 197;
% $$$ p  = 1;
% $$$ s1 = 4;
% $$$ 
% $$$ pfs =     MTAAknnpfs(Trial,units,states{s1},0,'numIter',1000, ...
% $$$                      'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
% $$$ 
% $$$ xs = sq(apfs.patchRateInd(u,s1,p,1,~isnan(sq(apfs.patchRateInd(u,s1,p,1,:)))&sq(apfs.patchRateInd(u,s1,p,1,:))~=0))
% $$$ ys = sq(apfs.patchRateInd(u,s1,p,2,~isnan(sq(apfs.patchRateInd(u,s1,p,2,:)))&sq(apfs.patchRateInd(u,s1,p,2,:))~=0))
% $$$ patchRatePos = [pfs.adata.bins{1}(xs),pfs.adata.bins{1}(ys)];
% $$$ 
% $$$ sigmap = pfs.plot(map(u,3),'sig');
% $$$ figure,pfs.plot(map(u,3));hold on,plot(apfs.patchCOM(u,s1,p,1),apfs.patchCOM(u,s1,p,2),'wo')
% $$$ plot(patchRatePos(:,1),patchRatePos(:,2),'w.')
% $$$ figure,imagesc(sigmap'<.05),axis xy
% $$$ 
% $$$ figure,imagesc(sigmap'),axis xy


npatch = 2;
minEl = 8;
niter = 10000;

% apa - all patches over lap

apa = repmat(struct('ovIntShStPerm',nan([numsts,numsts,npatch,npatch,niter]),'ovUniShStPerm',nan([numsts,numsts,npatch,niter])),size(map,1),1);

%apa.ovUniShStPerm = nan([numel(units),numsts,numsts,npatch,npatch,niter]);
%apa.ovIntShStPerm = nan([numsts,numsts,npatch,npatch,niter]);
%apa.pvpallShStPerm = nan([numel(units),numsts,numsts,npatch,npatch,niter]);

for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
    Trial.stc.updateMode(stc_mode);Trial.stc.load;
    for i = 1:numsts,
        pfs{ses,i} =     MTAAknnpfs(Trial,map(map(:,2)==ses,3),states{i},0,'numIter',1, ...
                            'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
    end
end

parfor u = map(:,1)',    
    for s1 = 1:numsts,
        for s2 = 1:numsts,
            for p1 = 1:npatch,
                for p2 = 1:npatch,

                    %------overlaped patch intersection shuffled state permutations------%
                    oLapInd1 = ismember(sq(apfs.patchRateInd(u,s1,p1,:,~isnan(sq(apfs.patchRateInd(u,s1,p1,1,:)))))',...
                                        sq(apfs.patchRateInd(u,s2,p2,:,~isnan(sq(apfs.patchRateInd(u,s2,p2,1,:)))))','rows');
                    if ~isempty(oLapInd1),
                        oLapInd2 = ismember(sq(apfs.patchRateInd(u,s2,p2,:,~isnan(sq(apfs.patchRateInd(u,s2,p2,1,:)))))',...
                                            sq(apfs.patchRateInd(u,s1,p1,:,~isnan(sq(apfs.patchRateInd(u,s1,p1,1,:)))))','rows');
                        if ~isempty(oLapInd2),
                            rm1 = sq(apfs.patchRateMap(u,s1,p1,~isnan(sq(apfs.patchRateMap(u,s1,p1,:)))));
                            rm2 = sq(apfs.patchRateMap(u,s2,p2,~isnan(sq(apfs.patchRateMap(u,s2,p2,:)))));
                            
                            if sum(oLapInd1)>(numel(oLapInd1)/2)&numel(oLapInd1)>minEl...
                                    & sum(oLapInd2)>(numel(oLapInd2)/2)&numel(oLapInd2)>minEl,
                                for i = 1:niter
                                    apa(u).ovIntShStPerm(s1,s2,p1,p2,i) = mean([randsample(rm1(oLapInd1),minEl/2,true);randsample(rm2(oLapInd2),minEl/2,true)]);
                                end
                            end
                        end
                    end

                    %------overlaped patch union shuffled state permutations------%
                    rmap_st1 = sq(apfs.patchRateMap(u,s1,p1,~isnan(sq(apfs.patchRateMap(u,s1,p1,:)))));
                    if ~isempty(rmap_st1),
                        rmap_st2 = pfs{map(u,2),s2}.plot(map(u,3));
                        patchPos_st1 = sq(apfs.patchRateInd(u,s1,p1,:,~isnan(sq(apfs.patchRateInd(u,s1,p1,1,:)))))';
                        if ~isempty(patchPos_st1),
                            rmap_st2 = rmap_st2(sub2ind(size(rmap_st2),patchPos_st1(:,1),patchPos_st1(:,2)));

                            for i = 1:niter
                                apa(u).ovUniShStPerm(s1,s2,p1,i) = mean([randsample(rmap_st1,minEl/2,true);randsample(rmap_st2,minEl/2,true)]);
                            end
                        end
                    end
                    %------end------%                    
                end
            end
        end
    end
end


parfor i = 1:1000,
    for j= 1:100,
    tstruct(i).tv(j)= i*j+i+j;

end
end


u = 216;
s1= 2;
s2= 3;
p1= 1;
p2= 1;

figure,hist(sq(apa(u,s1,s2,p1,p2,:)),30)

sq(apfs.patchRateMap(u,s1,p1,~isnan(sq(apfs.patchRateMap(u,s1,p1,:)))))

figure,

plot(1./(sum(sq(sort(apa(~isnan(apa(:,2,3,1,1,1)),2,3,1,1,:),2))<repmat(apfs.patchMFR(~isnan(apa(:,2,3,1,1,1)),3,1),[1,10000]),2)),1./(sum(sq(sort(apa(~isnan(apa(:,2,3,1,1,1)),2,3,1,1,:),2))<repmat(apfs.patchMFR(~isnan(apa(:,2,3,1,1,1)),2,1),[1,10000]),2)),'.');



s1= 4;
s2= 5;
figure,
for p1=1:2;
for p2=1:2;

nind = ~isnan(apa(:,s1,s2,p1,p2,1));
sig1 = 0.05<1./(sum(sq(sort(apa(nind,s1,s2,p1,p2,:),2))<repmat(apfs.patchMFR(nind,s2,1),[1,10000]),2));
sig2 = 0.05>1./(sum(sq(sort(apa(nind,s1,s2,p1,p2,:),2))<repmat(apfs.patchMFR(nind,s1,1),[1,10000]),2));

nmfr1 = apfs.patchMFR(nind,s1,1);
nmfr2 = apfs.patchMFR(nind,s2,1);
%nmfr1 = apfs.patchPFR(nind,s1,1);
%nmfr2 = apfs.patchPFR(nind,s2,1);
plot(nmfr1(sig1&sig2),nmfr2(sig1&sig2) ,'.c');
hold on;
plot(nmfr1(~sig1&~sig2),nmfr2(~sig1&~sig2) ,'.m');
end
end

%plot(nmfr1(~sig1&sig2),nmfr2(~sig1&sig2) ,'.g');
%hold on;
%plot(nmfr1(sig1&~sig2),nmfr2(sig1&~sig2) ,'.r');


sesId = 1
% $$$     if map(u,2)~=sesId,
% $$$         sesId = map(u,2);
% $$$         pfs = {};
% $$$         Trial = MTATrial(sesList{sesId}{1},sesList{sesId}{3},sesList{sesId}{2});
% $$$         Trial.stc.updateMode(stc_mode);Trial.stc.load;
% $$$         for i=1:numsts,
% $$$             pfs{i} = MTAAknnpfs(Trial,units,states{i},0,'numIter',1, ...
% $$$                                 'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
% $$$         end
% $$$     end
