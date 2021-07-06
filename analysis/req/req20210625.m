%req20210625
%    Tags: decoding 
%    Status: Active
%    Type: Analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: General
%    Description: decoding with ratemaps and second order unit coactivation maps
%
% GIVEN multiple place fields
% n = number of cells
% m = max number of fields
% fieldPos is nan where m>1 when only one field exists
% fieldPos[n,m,2]

% fieldAct is a logical matrix representing currently active field
% fieldAct[n,m]

% currentPos is the current estimated position

% prob x given cell 1 and 2 are active

% for each spike of cell 1 search for spike of cell 2 with a specified time interval
%     for each of such instance found assign a spike to psudo cell(1,2) at the midpoint
%         between spike time of cell 1 and spike time of cell 2.
%
%
% outer sum? of ufr

Trial = Trials{7};
tunits = units{7};
tstc = stc{7};
spk = create( copy(Trial.spk), Trial, 1, '', tunits, 'deburst');

xyz = preproc_xyz(Trial,'trb');
resample(xyz,16);
sts = [tstc{'t-m-s',16}];
sts.cast('TimeSeries',xyz);
sts = sts.data;
hxyz = xyz(:,'nose',[1,2]);


% PLOTPF args -----------------
bound_lims = [-500,500;
              -500,500];
binDims = [20,20];
SmoothingWeights = [2.5,2.5];
type = 'xy';
posSampleRate = xyz.sampleRate;
%------------------------------

cocluStart = spk.map(end,1)+1;
coclu = [];
cores = [];
ratemap = zeros([50*50,numel(tunits),numel(tunits)]);
for u1 = 1:numel(tunits)
    for u2 = u1+1:numel(tunits)
        disp([tunits(u1),tunits(u2)]);
tic
        res1 = spk(tunits(u1));
        res2 = spk(tunits(u2));
        resn = zeros([numel(res1)*numel(res2),1]);
        k = 1;

        for r1 = 1:numel(res1)
            for r2 = 1:numel(res2)
                rdiff = res1(r1)-res2(r2);
                resn(k) = (double(sign(rdiff) ==-1)*(res1(r1)+abs(rdiff)/2) + ...
                           double(sign(rdiff) == 1)*(res2(r2)+abs(rdiff)/2) + ...
                           double(sign(rdiff) == 0)*(res2(r2))).*double(abs(rdiff)<0.05);
                k = k+1;
            end
        end

        resn = nonzeros(resn);

        if numel(resn)<10,
            cocluStart = cocluStart+1;
            continue;
        end
        
        resnx = round(resn*xyz.sampleRate)+1;
        resnx(resnx<1|resnx>numel(sts)) =[];
        resnx = resnx(sts(resnx));

        pos = [xyz(sts,'hcom',1),xyz(sts,'hcom',2)];
        spkpos = [xyz(resnx,'hcom',1),xyz(resnx,'hcom',2)];

        [ratemap(:,u1,u2),Bins] = PlotPF(Trial,spkpos,pos,binDims,SmoothingWeights,type,bound_lims,posSampleRate);

        cores = cat(1,cores,resn);
        coclu = cat(1,coclu,cocluStart*ones([numel(resn),1]));
        cocluStart = cocluStart+1;
toc        
    end
end


ratemaps = [];
for u1 = 1:numel(tunits)
    for u2 = u1+1:numel(tunits)
        ratemaps = cat(2,ratemaps,ratemap(:,u1,u2));
    end
end
ratemaps(~maskcirc(:)) = nan;

figure();
for u = find(sum(ratemaps,'omitnan')>100),
    clf();
pax = pcolor(dc2.pfs{1}.adata.bins{:},(reshape(ratemaps(:,u),50,50).*maskcirc)');
pax.EdgeColor = 'none';
axis('xy');
colorbar();

cluOffset = 184;

resnx = round(cores(coclu==(cluOffset+u))*xyz.sampleRate)+1;
resnx(resnx<1|resnx>numel(sts)) =[];
resnx = resnx(sts(resnx));

hold('on')
lax = plot(xyz(resnx,'hcom',1),xyz(resnx,'hcom',2),'.r');;
xlim([-500,500]);
waitforbuttonpress()
end
cocluIds = 185:185+size(ratemaps,2);


selectedCocluIds = find(max(ratemaps)>1.5)+184;
%selectedCocluIds = find(sum(ratemaps,'omitnan')>300)+184;
selectedRatemaps = ratemaps(:,selectedCocluIds-184);
selectedCoclu = coclu(ismember(coclu,selectedCocluIds));
selectedCores = cores(ismember(coclu,selectedCocluIds));

pfs = dc2.pfs{7};

sts = [tstc{'t-m-s',16}];
sts.cast('TimeSeries',xyz);
sts = sts.data;

overwrite = true;
spkWindow = 0.3;
tag = ['decode_1Order_xy',...
       '_sr_',num2str(xyz.sampleRate),...
       '_sw_',num2str(spkWindow)];
cospk = copy(spk);
cospk.map = [];
cospk.clu = [spk.clu];
cospk.res = [spk.res];
allClu = [tunits];
cospk.res = round(cospk.res*xyz.sampleRate)+1;
cospk.clu(cospk.res<1|cospk.res>numel(sts)) =[];
cospk.res(cospk.res<1|cospk.res>numel(sts)) =[];
cospk.sampleRate = xyz.sampleRate;
coufr = Trial.load('ufr',xyz,cospk,allClu,spkWindow,'gauss',true);
copfs.adata.bins = Bins;
copfs.data.rateMap = pfs.data.rateMap;    
[cocom, comax, cosax, copost] =                                  ...
    decode_ufr(Trial,                                            ...
               allClu,                                           ...
               cospk.sampleRate,                                 ...
               coufr, copfs, [],                                 ...
               dc2.mask,                                         ...
               dc2.smoothingWeights,                             ...
               'tag',tag,                                        ...
               'overwrite',overwrite);

tag = ['decode_2Order_xy',...
       '_sr_',num2str(xyz.sampleRate),...
       '_sw_',num2str(spkWindow)];
caspk = copy(spk);
caspk.map = [];
caspk.clu = [spk.clu;selectedCoclu];
caspk.res = [spk.res;selectedCores];
allClu = [tunits,selectedCocluIds];
caspk.res = round(caspk.res*xyz.sampleRate)+1;
caspk.clu(caspk.res<1|caspk.res>numel(sts)) =[];
caspk.res(caspk.res<1|caspk.res>numel(sts)) =[];
caspk.sampleRate = xyz.sampleRate;
caufr = Trial.load('ufr',xyz,caspk,allClu,spkWindow,'gauss',true);
capfs.adata.bins = Bins;
capfs.data.rateMap = cat(2,pfs.data.rateMap,selectedRatemaps);
[cacam, camax, casax, capost] =                                  ...
    decode_ufr(Trial,                                            ...
               allClu,                                           ...
               caspk.sampleRate,                                 ...
               caufr, capfs, [],                                 ...
               dc2.mask,                                         ...
               dc2.smoothingWeights,                             ...
               'tag',tag,                                        ...
               'overwrite',overwrite);

sts = [tstc{'t-m-s-r',16}];
%sts = [tstc{'x+p&t',16}];
sts.cast('TimeSeries',xyz);
sts = sts.data;


ind = find((capost>0.06&copost>0.02)       ...
           & sts                            ...
           & sum(double(caufr.data>0.5),2)>3);
%ind = find((capost>0.01&copost>0.01)&sts&sum(double(coufr.data>0.5),2)>3);
nxyz = sq(xyz(:,'nose',1:2));

ind = ind(1:2:end);
%ind = ind(randn(size(ind))<0);

figure();
hold('on');
histogram(log10(sqrt(sum((cosax(ind,:)-nxyz(ind,:)).^2,2))./10),linspace(-1,2,50));
Lines(median(log10(sqrt(sum((cosax(ind,:)-nxyz(ind,:)).^2,2))./10)),[],'r');
Lines(mean(log10(sqrt(sum((cosax(ind,:)-nxyz(ind,:)).^2,2))./10)),[],'k');
mean(log10(sqrt(sum((cosax(ind,:)-nxyz(ind,:)).^2,2))./10))

histogram(log10(sqrt(sum((casax(ind,:)-nxyz(ind,:)).^2,2))./10),linspace(-1,2,50));
Lines(median(log10(sqrt(sum((casax(ind,:)-nxyz(ind,:)).^2,2))./10)),[],'r');
Lines(mean(log10(sqrt(sum((casax(ind,:)-nxyz(ind,:)).^2,2))./10)),[],'k');
mean(log10(sqrt(sum((casax(ind,:)-nxyz(ind,:)).^2,2))./10))
grid('on');



