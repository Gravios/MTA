%function dstruct = population_placeFieldAnalysis(sesList,states,pftype)
MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute'); 

sesList = {{'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};

states = {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'};

pftype = 'MTAAknnpf';

numsts = numel(states);
pfstats = {};
    pfs={};
for ses = 2:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
    Trial.xyz.load(Trial);
    Trial.load('nq');

    %units = find(Trial.nq.SpkWidthR>0.7&Trial.nq.eDist>30);
    units = [];

    % Load PlaceFields

    for i = 1:numsts,
        switch pftype
          case 'MTAAknnpf'
            try
            pfs{ses,i} = MTAAknnpfs(Trial,units,states{i},0,'numIter',1000, ...
                            'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
            end
          case 'MTAApfs'
            % not ready 
            % pfs{i} = MTAApfs(Trial,units,states{i},0,'numIter',1000)
        end

        if i==1,
            units = pfs{ses,1}.data.clu;
        end
    end



    % Calulate the characteristic values of place fields
    for u = units
        for state = 1:numsts
            pfstats{ses}(u==units,state) = PlaceFieldStats(Trial,pfs{ses,state},u);    
        end
    end


tic
    % Initialize correlation vars
    stsCor = zeros(numel(units),numsts,numsts,2,2,2);
    % Initialize overlap vars
    overlap = zeros(numel(units),numsts,numsts);
    overlapMFR = zeros(numel(units),numsts,numsts,2);
    % Initialize Permutation test vars
    numIter=1000;
    permsamp = zeros(numel(units),numsts,numsts,numIter);
toc
tic
    for u = units,
        for statei = 1:numsts,
            for statej = 1:numsts,

                hsig = pfs{ses,statei}.plot(u,'sig');
                lsig = pfs{ses,statej}.plot(u,'sig');
                hmap = pfs{ses,statei}.plot(u);
                lmap = pfs{ses,statej}.plot(u);
                mind = hsig<.05&lsig<.05;

                if sum(mind(:))>10, 
                    % state ratemap correlations
                    [stsCor(u==units,statei,statej,:,:,1),...
                     stsCor(u==units,statei,statej,:,:,2)] = corrcoef(hmap(mind),lmap(mind));
                    % States ratemap overlap
                    novlp = sum(mind(:))/(sum(reshape(hsig<.05,[],1))+sum(reshape(hsig<.05,[],1)));
                    overlap(u==units,statei,statej) = novlp;
                    overlapMFR(u==units,statei,statej,:) = [mean(hmap(mind)),mean(lmap(mind))];
                    perminds = randi(sum(mind(:))*2,[sum(mind(:)),numIter]);
                    permpool = [hmap(mind);lmap(mind)];
                    permsamp(u==units,statei,statej,:) = sq(mean(permpool(perminds)));
                else 
                    stsCor(u==units,statei,statej,:,:,2) = ones(2,2);
                end
                
            end

        end
    end
toc
    dstruct(ses).clu   = pfs{ses,1}.data.clu(:);
    dstruct(ses).elClu = pfs{ses,1}.data.elClu(:);
    dstruct(ses).el    = pfs{ses,1}.data.el(:);
    dstruct(ses).ratecorr = stsCor;   
    dstruct(ses).permsamp = permsamp; 
    dstruct(ses).overlap  = overlap;
    dstruct(ses).overlapMFR  = overlapMFR;
    nqfields = fieldnames(Trial.nq);
    for i = 1:numel(nqfields),
        dstruct(ses).(nqfields{i}) = Trial.nq.(nqfields{i})(units,:);  
    end
end

cstruct = CatStruct(dstruct,[],1);

cpfstats = cat(1,pfstats{:});


%% Distance vs FR difference
s1 = 2;
s2 = 3;
%cind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>20&cstruct.FirRate>0.2&[[cpfstats(:,s1).patchArea]>2000|[cpfstats(:,s2).patchArea]>2000]');
cind = find(cstruct.SpkWidthR>0.3&cstruct.TimeSym>2.5&cstruct.eDist>19&cstruct.FirRate>0.2);

dist = sqrt(sum((reshape([cpfstats(cind,s1).patchCOM],2,[])'-reshape([cpfstats(cind,s2).patchCOM],2,[])').^2,2));
rdr = ([cpfstats(cind,s1).peakFR]-[cpfstats(cind,s2).peakFR])./ ...
      ([cpfstats(cind,s1).peakFR]+[cpfstats(cind,s2).peakFR]);
mrdr = ([cpfstats(cind,s1).patchMFR]-[cpfstats(cind,s2).patchMFR])./ ...
      ([cpfstats(cind,s1).patchMFR]+[cpfstats(cind,s2).patchMFR]);

pr = max([[cpfstats(cind,s1).peakFR];[cpfstats(cind,s2).peakFR]])';

pra = max([[cpfstats(:,s1).peakFR];[cpfstats(:,s2).peakFR]])';

figure,plot(dist,mrdr,'d','MarkerFaceColor','b')

figure,scatter(dist,mrdr,pr.*5);

figure,hist(mrdr,30,'MarkerFaceColor','b')
figure,hist(mrdr(pr>4&pr<40),30,'MarkerFaceColor','b')

figure,hist(cstruct.ratecorr(cstruct.overlap(:,2,3)>0.05&pra>4&pra<40&cstruct.SpkWidthR>0.3&cstruct.TimeSym>2.5&cstruct.eDist>19&cstruct.FirRate>0.2,2,3,1,2,1),15)

figure,plot([cpfstats(cind,s1).peakFR]-[cpfstats(cind,s2).peakFR],[cpfstats(cind,s1).peakFR]+[cpfstats(cind,s2).peakFR],'d','MarkerFaceColor','b')
line([0;-30],[0;30])
line([0;30],[0;30])
line([0;30],[0;90],'color',[0,0,1])
line([0;-30],[0;90],'color',[0,0,1])
ylim([0,80])
line([0;0],[0;80],'color',[0,0,0])
hold on,plot([cpfstats(rind,s1).peakFR]-[cpfstats(rind,s2).peakFR],[cpfstats(rind,s1).peakFR]+[cpfstats(rind,s2).peakFR],'rd','MarkerFaceColor','r')


figure,plot([cpfstats(cind,s1).patchMFR]-[cpfstats(cind,s2).patchMFR],[cpfstats(cind,s1).patchMFR]+[cpfstats(cind,s2).patchMFR],'d','MarkerFaceColor','b')
line([0;20],[0;20])
line([0;-20],[0;20])
line([0;20],[0;60],'color',[0,0,1])
line([0;-20],[0;60],'color',[0,0,1])
ylim([0,60])
line([0;0],[0;60],'color',[0,0,0])
hold on,plot([cpfstats(rind,s1).patchMFR]-[cpfstats(rind,s2).patchMFR],[cpfstats(rind,s1).patchMFR]+[cpfstats(rind,s2).patchMFR],'rd','MarkerFaceColor','r')


% $$$ line([0;-20],[0;60],'color',[1,0,1])
% $$$ line([0;20],[0;60],'color',[1,0,1])
% $$$ 
% $$$ line([0;20],[0;80],'color',[1,0,0])
% $$$ line([0;-20],[0;80],'color',[1,0,0])
% $$$ 
% $$$ line([0;20],[0;100],'color',[1,0,0])
% $$$ line([0;-20],[0;100],'color',[1,0,0])
% $$$ 



s1 = 2;
s2 = 3;
rslope = ([cpfstats(:,s1).patchMFR]+[cpfstats(:,s2).patchMFR])./([cpfstats(:,s1).patchMFR]-[cpfstats(:,s2).patchMFR]);
rslope = ([cpfstats(:,s1).peakFR]+[cpfstats(:,s2).peakFR])./([cpfstats(:,s1).peakFR]-[cpfstats(:,s2).peakFR]);
minrate = [cpfstats(:,s1).peakFR]>5|[cpfstats(:,s2).peakFR]>5;
maxrate = [cpfstats(:,s1).peakFR]<30|[cpfstats(:,s2).peakFR]<40;
 
rind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>25&cstruct.FirRate>0.2&([rslope<3&rslope>0&minrate&maxrate]'));
rind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>25&cstruct.FirRate>0.2&([rslope>-3&rslope<0&minrate&maxrate]'));
rind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>25&cstruct.FirRate>0.2&([rslope<-3&minrate&maxrate]'));
rind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>25&cstruct.FirRate>0.2&([rslope>3&minrate&maxrate]'));

rind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>25&cstruct.FirRate>0.2&([(rslope<-3|rslope>3)&minrate&maxrate]'));
for i = rind(:)'
clf
pfr = pfs{cstruct.ses(i),s1}.plot(cstruct.clu(i));
pfw = pfs{cstruct.ses(i),s2}.plot(cstruct.clu(i));
clim =max([max(pfr),max(pfw)]);
name = [strjoin(sesList{cstruct.ses(i)},'.'), '-unit-' num2str(cstruct.clu(i))];
set(gcf,'Name',name)
subplot(211);hold on
pfs{cstruct.ses(i),s1}.plot(cstruct.clu(i),[],[],[0,clim]);
%set(gca,'XTickLabel',{}),set(gca,'YTickLabel',{})
title(states{s1})
%circle(0,0,410,'m'),
subplot(212);hold on
pfs{cstruct.ses(i),s2}.plot(cstruct.clu(i),[],[],[0,clim]);
%set(gca,'XTickLabel',{}),set(gca,'YTickLabel',{})
title(states{s2})
%circle(0,0,410,'m'),
print(gcf,'-depsc',['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/similar_rear_walk/' name '.eps'])
print(gcf,'-dpng',['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/similar_rear_walk/' name '.png'])
end

rdist = sqrt(sum((reshape([cpfstats(rind,s1).patchCOM],2,[])'-reshape([cpfstats(rind,s2).patchCOM],2,[])').^2,2));
figure,plot(rdist,cstruct.ratecorr(rind,2,3,1,2,1),'d');

s1 = 4;
s2 = 5;
wslope = ([cpfstats(:,s1).peakFR]+[cpfstats(:,s2).peakFR])./([cpfstats(:,s1).peakFR]-[cpfstats(:,s2).peakFR]);
wind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>25&cstruct.FirRate>0.2&([wslope>-3&wslope<0&minrate&maxrate]'));
wind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>25&cstruct.FirRate>0.2&([wslope<3&wslope>0&minrate&maxrate]'));

wind = find(cstruct.SpkWidthR>0.3&cstruct.eDist>25&cstruct.FirRate>0.2&([wslope>3&minrate&maxrate]'));
for i = wind(:)'
clf
pfr = pfs{cstruct.ses(i),s1}.plot(cstruct.clu(i));
pfw = pfs{cstruct.ses(i),s2}.plot(cstruct.clu(i));
clim =max([max(pfr),max(pfw)]);
name = [strjoin(sesList{cstruct.ses(i)},'.'), '-unit-' num2str(cstruct.clu(i))];
set(gcf,'Name',name)
subplot(211);hold on
pfs{cstruct.ses(i),s1}.plot(cstruct.clu(i),[],[],[0,clim]);
%set(gca,'XTickLabel',{}),set(gca,'YTickLabel',{})
title(states{s1})
%circle(0,0,410,'m'),
subplot(212);hold on
pfs{cstruct.ses(i),s2}.plot(cstruct.clu(i),[],[],[0,clim]);
%set(gca,'XTickLabel',{}),set(gca,'YTickLabel',{})
title(states{s2})
%circle(0,0,410,'m'),
print(gcf,'-depsc',['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/lwalk_hwalk/' name '.eps'])
print(gcf,'-dpng',['/gpfs01/sirota/bach/homes/gravio/figures/population_pfs/lwalk_hwalk/' name '.png'])
end
