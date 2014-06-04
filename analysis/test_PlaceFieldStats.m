
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

apa = repmat(struct('ovIntShStPerm',nan([numsts,numsts,npatch,npatch,niter]),...
                    'ovUniShStPerm',nan([numsts,numsts,npatch,niter]),...
                    'ovUniShStPri', nan([numsts,numsts,npatch,niter]),...
                    'ovUniShStSec', nan([numsts,numsts,npatch,niter])...
             ),size(map,1),1);

%apa.ovUniShStPerm = nan([numel(units),numsts,numsts,npatch,npatch,niter]);
%apa.ovIntShStPerm = nan([numsts,numsts,npatch,npatch,niter]);
%apa.pvpallShStPerm = nan([numel(units),numsts,numsts,npatch,npatch,niter]);

pfs = {};
for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
    Trial.stc.updateMode(stc_mode);Trial.stc.load;
    for i = 1:numsts,
        pfs{ses,i} =     MTAAknnpfs(Trial,map(map(:,2)==ses,3),states{i},0,'numIter',1, ...
                            'ufrShufBlockSize',0,'binDims',[20,20],'distThreshold',70);
    end
end

parfor u = map(1:221,1)',    
    for s1 = 1:numsts,
        for s2 = 1:numsts,
            for p1 = 1:npatch,
                for p2 = 1:npatch,

                    %------overlaped patch intersection shuffled state permutations------%
% $$$                     oLapInd1 = ismember(sq(apfs.patchRateInd(u,s1,p1,:,~isnan(sq(apfs.patchRateInd(u,s1,p1,1,:)))))',...
% $$$                                         sq(apfs.patchRateInd(u,s2,p2,:,~isnan(sq(apfs.patchRateInd(u,s2,p2,1,:)))))','rows');
% $$$                     if ~isempty(oLapInd1),
% $$$                         oLapInd2 = ismember(sq(apfs.patchRateInd(u,s2,p2,:,~isnan(sq(apfs.patchRateInd(u,s2,p2,1,:)))))',...
% $$$                                             sq(apfs.patchRateInd(u,s1,p1,:,~isnan(sq(apfs.patchRateInd(u,s1,p1,1,:)))))','rows');
% $$$                         if ~isempty(oLapInd2),
% $$$                             rm1 = sq(apfs.patchRateMap(u,s1,p1,~isnan(sq(apfs.patchRateMap(u,s1,p1,:)))));
% $$$                             rm2 = sq(apfs.patchRateMap(u,s2,p2,~isnan(sq(apfs.patchRateMap(u,s2,p2,:)))));
% $$$                             
% $$$                             if sum(oLapInd1)>(numel(oLapInd1)/2)&numel(oLapInd1)>minEl...
% $$$                                     & sum(oLapInd2)>(numel(oLapInd2)/2)&numel(oLapInd2)>minEl,
% $$$                                 for i = 1:niter
% $$$                                     apa(u).ovIntShStPerm(s1,s2,p1,p2,i) = mean([randsample(rm1(oLapInd1),minEl/2,true);randsample(rm2(oLapInd2),minEl/2,true)]);
% $$$                                 end
% $$$                             end
% $$$                         end
% $$$                     end

                    %------overlaped patch union shuffled state permutations------%
                    ind = sq(apfs.patchRateMap(u,s1,p1,:))~=0&~isnan(sq(apfs.patchRateMap(u,s1,p1,:)));
                    rmap_st1 = sq(apfs.patchRateMap(u,s1,p1,ind));
                    if ~isempty(rmap_st1),
                        rmap_st2 = pfs{map(u,2),s2}.plot(map(u,3));
                        patchPos_st1 = sq(apfs.patchRateInd(u,s1,p1,:,ind))';
                        if ~isempty(patchPos_st1),
                            rmap_st2 = rmap_st2(sub2ind(size(rmap_st2),patchPos_st1(:,1),patchPos_st1(:,2)));
                            if numel(rmap_st1)>=minEl&numel(rmap_st2)>=minEl
                            for i = 1:niter
                                apa(u).ovUniShStPri(s1,s2,p1,i) = mean([randsample(rmap_st1,minEl/2,true);randsample(rmap_st1,minEl/2,true)]);
                                apa(u).ovUniShStSec(s1,s2,p1,i) = mean([randsample(rmap_st2,minEl/2,true);randsample(rmap_st2,minEl/2,true)]);
                                apa(u).ovUniShStPerm(s1,s2,p1,i) = mean([randsample(rmap_st1,minEl/2,true);randsample(rmap_st2,minEl/2,true)]);
                            end
                            end
                        end
                    end
                    %------end------%                    
                end
            end
        end
    end
end

ous.perm = permute(cat(5,apa.ovUniShStPerm),[4,5,1,2,3]);
ous.pri = permute(cat(5,apa.ovUniShStPri),[4,5,1,2,3]);
ous.sec = permute(cat(5,apa.ovUniShStSec),[4,5,1,2,3]);

parfor u = map(1:221,1)',    
    for st1 = 1:numsts,
        for st2 = 1:numsts,
            for p1 = 1:npatch,
                    pmean = mean(ous.perm(:,u,st1,st2,p1));
                    pstd  =std(ous.perm(:,u,st1,st2,p1));
                    [pri(u).h(st1,st2,p1),pri(u).p(st1,st2,p1),pri(u).ci(st1,st2,p1,:),pri(u).zval(st1,st2,p1)] ...
                        =ztest(ous.pri(:,u,st1,st2,p1),pmean,pstd);
                    [sec(u).h(st1,st2,p1),sec(u).p(st1,st2,p1),sec(u).ci(st1,st2,p1,:),sec(u).zval(st1,st2,p1)] ...
                        =ztest(ous.sec(:,u,st1,st2,p1),pmean,pstd);

            end
        end
    end
end

sec_z = permute(cat(5,sec.zval),[5,1,2,3,4]);
pri_z = permute(cat(5,pri.zval),[5,1,2,3,4]);





ous.perm = permute(cat(5,apa.ovUniShStPerm),[4,5,1,2,3]);
ous.pri = permute(cat(5,apa.ovUniShStPri),[4,5,1,2,3]);
ous.sec = permute(cat(5,apa.ovUniShStSec),[4,5,1,2,3]);



figure,plot(pri_z(:,2,5,1),pri_z(:,2,4,1),'.')

line([-6000;6000],[-6000;6000])


s1=2;s2=4;

%state patch distance
spd = sqrt(sum(sq((apfs.patchCOM(:,s1,1,:))-apfs.patchCOM(:,s2,1,:)).^2,2));
%relative rate difference
rrd = sq((apfs.patchMFR(:,s1,1)-apfs.patchMFR(:,s2,1))./(apfs.patchMFR(:,s1,1)+apfs.patchMFR(:,s2,1)));
%rrd = sq((apfs.patchPFR(:,s1,1)-apfs.patchPFR(:,s2,1))./(apfs.patchPFR(:,s1,1)+apfs.patchPFR(:,s2,1)));
%rrd = sq((apfs.peakFR(:,s1,1)-apfs.peakFR(:,s2,1))./(apfs.peakFR(:,s1,1)+apfs.peakFR(:,s2,1)));

figure,hist(rrd(spd<100),20);
figure,hist(rrd(:),20);

figure,plot(spd,rrd,'.')
