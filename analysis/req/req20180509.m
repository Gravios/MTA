% Tasks 20180509
%
% 1. Compile: rate, size, position, and spatial information, of each place field state
% 2. plot: state x state patch distance 


% Compute shuffled state pfs
MjgER2016_load_data;



tind = 20; % jg05-20120312.cof.all


pfs = pfs_2d_states(Trials{tind},units{tind},[],[],false,'',false);
pfsArgs = struct('posShuffle',true);
pfsShuffled = pfs_2d_states(Trials{tind},units{tind},[],[],false,'',false,pfsArgs);

pft = pfs_2d_theta(Trials{tind},units{tind});
pftArgs = struct('posShuffle',true);
pftShuffled = pfs_2d_theta(Trials{tind},units{tind},false,false,pftArgs);


pfs = cat(2,{pft},pfs);
pfsShuffled = cat(2,{pftShuffled},pfsShuffled);

numStates = numel(pfs);

%pbins = interpParPfs.bins; ipp = interpParPfs;
pbins = pfs{1}.adata.bins; ipp = [];


bins = cell([1,2]);
[bins{:}] = ndgrid(pbins{:});

for u = 1:numel(units{tind});
clf();
mrate = max(cell2mat(cf(@(p,u) max(max(plot(p,u,'mean',false,[],true,0.5,false,[]))), pfs,repmat({units{tind}(u)},size(pfs)))));
for s = 1:numStates;

    rmapC = plot(pfs{s},units{tind}(u),'mean',false,[],true,0.5,false,ipp);
    rmapSm = plot(pfsShuffled{s},units{tind}(u),'mean',false,[],true,0.5,false,ipp);
    rmapSs = plot(pfsShuffled{s},units{tind}(u),'std',false,[],true,0.5,false,ipp);
% $$$     rmapSm = plot(pftShuffled,units{tind}(u),'mean',false,[],true,0.5,false,ipp);
% $$$     rmapSs = plot(pftShuffled,units{tind}(u),'std',false,[],true,0.5,false,ipp);
% $$$     rmapSm = plot(pfs{s},units{tind}(u),'mean',false,[],true,0.5,false,ipp);
% $$$     rmapSs = plot(pfs{s},units{tind}(u),'std',false,[],true,0.5,false,ipp);
    rmapSnr = plot(pfs{s},units{tind}(u),'snr',false,[],true,0.5,false,ipp);
    

    subplot2(numStates,5,s,1);  imagescnan({pbins{:},rmapC'},[0,mrate],[],true,[],[],[],@jet); 
    hold('on');                 contour(bins{:},rmapC,repmat(max(rmapC(:)),[1,2]).*0.5,'m')
                                title(num2str(units{tind}(u)));
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});
                            
    subplot2(numStates,5,s,2);  imagescnan({pbins{:},rmapSm'},[],[],true,[],[],[],@jet); 
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});    
    subplot2(numStates,5,s,3);  imagescnan({pbins{:},rmapSs'},[],[],true,[],[],[],@jet); 
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});


    zmap = (rmapC'-rmapSm')./rmapSs';
    subplot2(numStates,5,s,4);  imagescnan({pbins{:},zmap},[0,10],[],true,[],[],[],@jet);
    hold('on');                 contour(bins{:},zmap',repmat(3,[1,2]),'m')
                                title(pfs{s}.parameters.states);
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});
                               
    subplot2(numStates,5,s,5);  imagescnan({pbins{:},rmapSnr'},[0,10],[],true,[],[],[],@jet); 
    hold('on');                 contour(bins{:},rmapSnr,repmat(3,[1,2]),'m')
                                title(pfs{s}.parameters.states);
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});
end
waitforbuttonpress();
end




% $$$ cf(@(t,u,a) pfs_2d_states(t,u,[],[],false,'',false,a), Trials,units,repmat({pfsArgs},size(Trials)));
% $$$ cf(@(t,u,a) pfs_2d_theta(t,u,false,false,a),           Trials,units,repmat({pfsArgs},size(Trials)));

pfstats  = {};
pfbstats = {};
pfmstats = {};
for tind = 1:numel(Trials),
pfs = pfs_2d_states(Trials{tind},units{tind},[],[],false,'',false);
pft = pfs_2d_theta(Trials{tind},units{tind});
pfs = cat(2,{pft},pfs);
pftArgs = struct('posShuffle',true);
pftShuffled = pfs_2d_theta(Trials{tind},units{tind},false,true,pftArgs);
for u = 1:numel(units{tind}),
    for s = 1:numel(pfs),
        [pfstats{tind}{u}{s},pfbstats{tind}{u}{s},pfmstats{tind}{u}{s}] = ...
            compute_placefield_stats(Trials{tind},pfs{s},units{tind}(u),2,true,pftShuffled);
    end    
end
end



save('/storage/gravio/data/project/general/analysis/req20180509.mat','pfstats','pfbstats','pfmstats');


temp = [pfstats{:}]';
temp = cf(@(t) [t{:}],temp);
temp = vertcat(temp{:});

placefieldStats = {};
for s = 1:numel(pfs),
    placefieldStats{s} = CatStruct(temp(:,s),[],1);
end




% NUMBER of patches
patchCounts = sum(~isnan(placefieldStats{1}.patchPFR),2);
ind = ':';
eds = [-0.5,0.5,1.5,2.5];
figure,bar(eds,histc(patchCounts(ind),eds),'histc');




temp = [pfbstats{:}]';
temp = cf(@(t) [t{:}],temp);
temp = vertcat(temp{:});

pfStats = {};
for s = 1:numel(pfs),
    pfStats{s} = CatStruct(temp(:,s),[],1);
end

i = 3;
j = 4;

mfrWalk = mean(pfStats{i}.patchMFR(:,:,1),2,'omitnan');
mfrRear = mean(pfStats{j}.patchMFR(:,:,1),2,'omitnan');
ppDist  = median(sqrt(sum((median(pfStats{i}.patchCOM(:,:,:,1),3,'omitnan')...
                           -median(pfStats{j}.patchCOM(:,:,:,1),3,'omitnan')).^2,3)),2,'omitnan');


ppDist  = median(sqrt(sum((median(pfStats{i}.patchCOM(:,:,:,1),3,'omitnan')...
                           -median(pfStats{j}.patchCOM(:,:,:,1),3,'omitnan')).^2,3)),2,'omitnan');



mfrWalk(isnan(mfrWalk)) = 0;
mfrRear(isnan(mfrRear)) = 0;
ind = ppDist<100;

figure();
plot(mfrWalk(ind),...
     mfrRear(ind),...
     '.');
daspect([1,1,1]);
line([0,10],[0,10]);



for s = 1:numel(pfs),
    temp = pfbstats(:,s);
    pfbs{s} = cat(1,temp{:});
end
pfbs = cat(2,pfbs{:});



figure();
u = 3;
for s = 1:numel(pfs),
    subplot(numel(pfs),1,s);
    hist(pfbstats{u,s}.patchMFR(:,:,1),100)        
    xlim([0,10]);
end
pfbs = {};
for s = 1:numel(pfs),
    temp = pfbstats(:,s);
    pfbs{s} = cat(1,temp{:});
end
pfbs = cat(2,pfbs{:});




pbmr = reshape(mean(vertcat(pfbs.patchMFR),2,'omitnan'),[size(pfbstats),2]);

figure,imagesc(pbmr(:,[2:7],1)')

figure();
for s = 1:numel(pfs),
    subplot(numel(pfs),1,s);
    hist(pbmr(:,s,1),30)        
    xlim([0,10]);
end


figure();
hist(sqrt(sum(sq(pfbstats{u,s}.patchMFR(:,:,1,:)-pfbstats{u,1}.patchMFR(:,:,1,:)).^2,2)),100)


% Compute the placefield patch distance betwee theta and 
figure();
hist(sqrt(sum(sq(pfbstats{u,s}.patchMFR(:,:,1,:)-pfbstats{u,1}.patchMFR(:,:,1,:)).^2,2)),100)

figure;
plot(placefieldStats{1}.patchArea(:,1)./(40^2),...
     placefieldStats{1}.patchMFR(:,1),'.')

s = 1
figure;
plot(placefieldStats{s}.patchArea(:,1)./(40^2),...
     placefieldStats{s}.patchMFR(:,1),'.')


sum(~isnan(placefieldStats{5}.patchMFR(:,1)))


XMin = -pi;
XMax =  pi;
mycdf = @(x) cdf('norm',x,0,1);
probspace = @(CDF, XMin, XMax, N) arrayfun( @(p) max(fsolve( @(x) CDF(x)-p,...
                                                  [-0.01, 0.01],...
                                                  optimoptions('fsolve','Algorithm','trust-region'))), linspace(0,1,N) );
edy = probspace(mycdf,0.02,0.98,ny);


figure,
for j = 1:8,
    k = 1;
    for fields = fieldnames(placefieldStats{5})',
        subplot2(8,8,j,k);
        nind = nniz(placefieldStats{j}.(fields{1})(:,1));
        if strcmp('patchCOM',fields{1}),
            value = sqrt(sum(placefieldStats{j}.(fields{1})(nind,1,:).^2,3));
        else
            value = placefieldStats{j}.(fields{1})(nind,1);
        end
        
        if sum(nind)~=0,
        xbins = linspace([prctile(value,[2,98]),30]);
        bar(xbins,histc(value,xbins),'histc');
        end
        if j==1 ,title(fields{1}),end
        if k==1 ,ylabel(pfs{j}.parameters.states),end        
        k = k+1;
    end
end







% Low res multi-dimensional placefields
% $$$ 
% $$$ pfsArgs = struct(                                                                    ...
% $$$     'states',                         'theta-groom-sit',                             ...
% $$$     'overwrite',                      0,                                             ...
% $$$     'tag',                            [],                                            ...
% $$$     'binDims',                        [100,100],                                     ...
% $$$     'SmoothingWeights',               [0.8,0.8],                                     ...
% $$$     'type',                           'xy',                                          ...
% $$$     'spkShuffle',                     0,                                             ...
% $$$     'posShuffle',                     0,                                             ...
% $$$     'numIter',                        1,                                             ...
% $$$     'xyzp',                           MTADxyz([]),                                   ...
% $$$     'boundaryLimits',                 [-500,500;-500,500;-2,2;-2,2],                 ...
% $$$     'bootstrap',                      0,                                             ...
% $$$     'halfsample',                     0,                                             ...
% $$$     'shuffleBlockSize',               1,                                             ...
% $$$     'trackingMarker',                 'nose',                                        ...
% $$$     'autoSaveFlag',                   true,                                          ...
% $$$     'spkMode',                        'deburst'                                      ...
% $$$     );
% $$$ pfsArgs = struct2varargin(pfsArgs);
% $$$ pfn = MTAApfs(Trials{20},units{20},pfsArgs{:});
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ for u = 1:numel(units{20}),
% $$$     plot(pfn,units{20}(u),1,true,[],false);
% $$$     title(num2str(units{20}(u)));
% $$$     waitforbuttonpress();
% $$$ end
