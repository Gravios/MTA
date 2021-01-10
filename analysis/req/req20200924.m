% req20200924
%    Tags: ego, placefield, shifts, spatial information
%    Status: active
%    Type: Analysis
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: MjgER2016: BhvPlaceCode
%    Description: spatial information of egoPlaceFields as a function of placefield center shifts

MjgER2016_load_data();
Trial = Trials{20};
rot = 0.17;
hbaCorrection = -0.2;
overwrite = true;

sampleRate = 250;
pfsState = 'theta-groom-sit-rear';
spkMode = 'deburst';
binPhzs = linspace(-pi,pi,6);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
hbaBinEdges = -1.5:0.6:1.5;

phzCorrection = pi/4;
thetaPhzChan = sessionList(20).thetaRefGeneral;

units = [ 20, 21, 25, 31, 35, 44, 52, 61, 72, 79, 80, 81, 85,...% jg05-20120312  CA1
             103,104,110,111,116,138,139,151];    %109,

xyz = preproc_xyz(Trial,'trb');
xyz.filter('ButFilter',3,30,'low');
xyz.resample(sampleRate);

spk = Trial.load('spk',sampleRate,'gper',units,'deburst');

pft = pfs_2d_theta(Trial,units,'pfsArgsOverride', struct('halfsample',false,'numIter',1));



%function [pfs] = req20200924(Trial,units,xyz,spk,pft,rot,hbaCorrection,overwrite)

sampleRate = xyz.sampleRate;
binPhzs = linspace(0,2*pi,6);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
pfs = cell([1,numel(binPhzc)]);
verbose = true;

if verbose,
    disp(['[status]        compute_egohba_ratemap: processing trial: ',Trial.filebase]);
end

if isempty(units),
    return;
end;% if

% COMPUTE head basis
hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,...
                 [cos(rot),-sin(rot);sin(rot),cos(rot)],...
                 [2,3],...
                 [1,2]);

% GET theta state behaviors, minus rear
thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);

pfTemp = Trial;

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = units;
pargs.tag          = 'egofield';
pargs.binDims      = [20, 20]; % X Y 
pargs.SmoothingWeights = [2, 2];
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-400,400;-400,400];
pargs.states       = '';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
electrode = 0;

% CHECK existence of pfs object
shifts = -100:10:100;

pfs = {};

for x = 1:numel(shifts),
    for y = 1:numel(shifts),
        pargs.tag = ['egofield_refine_SW',num2str(pargs.SmoothingWeights(1)),...
                     '_X',num2str(shifts(x)),'_Y',num2str(shifts(y))];
        filepath = fullfile(Trial.spath, [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
        
        if exist(filepath,'file'),
            pfs{x,y} = load(filepath,'Pfs').('Pfs');
            if overwrite,
                pfs.purge_savefile();
            end;% if overwrite
        end;% if exist
        
        for unit = 1:numel(units),
            if unit==1 | electrode ~= spk.map(spk.map(:,1)==units(unit),2), % update phase state
                pargs.spk = copy( spk );
                electrode = 1;
                pargs.states = copy( thetaState );
                pargs.states.label = ['theta'];
                cast( pargs.states, 'TimePeriods' );
                resInd = WithinRanges( pargs.spk.res, pargs.states.data );
                pargs.spk.res = pargs.spk.res( resInd );
                pargs.spk.clu = pargs.spk.clu( resInd );
            end;%if unit==1
            
            [mxr,mxp] = pft.maxRate(units(unit));
            pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                              [multiprod(bsxfun(@minus,                        ...
                                                              mxp+shifts([x,y]),                ...
                                                              sq(xyz(:,'hcom',[1,2]))),      ...
                                                         hvec,2,[2,3])],           ...
                                               sampleRate,                                      ...
                                              'egocentric_placefield',                         ...
                                              'egopfs',                                        ...
                                              'p'                                              ...
                                              );
            
            pargs.xyzp = pfsCenterHR;
            pargs.units  = units(unit);
            pfsArgs = struct2varargin(pargs);
            pfTemp = MTAApfs(pfTemp,pfsArgs{:});
            if unit==1,
                try,
                    pfTemp.purge_savefile();
                end;
                pfTemp.save();        
            end;%if unit==1
        end;%for unit
        pfTemp.save();
        pfs{x,y} = pfTemp;    
        pfTemp = Trial;
        
    end;%for y
end;%for x


%%%<<< DISPLAY ratemaps with SI map
pfsa = decapsulate_and_concatenate_mtaapfs(pfs(:)',repmat({units},[1,prod(size(pfs))]));
pfsa = reshape(pfsa,[size(pfsa,1),numel(units),size(pfs,1),size(pfs,2)]);

xyi = [21,1;21,11;21,21;11,1;11,11;11,21;1,1;1,11;1,21];
figure();
for i= 1:9,
subplot(3,3,i);
pcolor( pfs{1}.adata.bins{:}, reshape(pfsa(:,3,xyi(i,1),xyi(i,2)),pfs{1}.adata.binSizes') );
caxis   ([0,8.2]);
colormap('jet');
shading ('flat');
axis    ('xy');
end

bins = pfs{1}.adata.bins;
binSizes = pfs{1}.adata.binSizes;
width = binSizes(1);
height =binSizes(2);
radius = round(binSizes(1)/2)-find(bins{1}<=-300,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mazeMask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);

pfsavb = pfsa(logical(mazeMask(:)),:,:,:);

psi = permute(sum(repmat(1./sum(~isnan(pfsavb(:,:,1))),[size(pfsavb,1),1]) ...
              .*bsxfun(@rdivide,pfsavb,mean(pfsavb,'omitnan')) ...
              .*log2(bsxfun(@rdivide,pfsavb,mean(pfsavb,'omitnan'))),'omitnan'),...
              [3,4,2,1]);

u = 1;
figure

mrate = 10;
figure();
for u = 1:21,
subplot(141);
plot(pft,units(u),1,'text',[0,mrate]);
subplot(142);
imagesc(shifts,shifts,psi(:,:,u)');
axis('xy');
subplot(143);
pcolor( pfs{1}.adata.bins{:}, reshape(pfsa(:,u,11,11),pfs{1}.adata.binSizes') );
caxis   ([0,mrate]);
colormap('jet');
shading ('flat');
axis    ('xy');
xlim([-300,300]);
ylim([-300,300]);
subplot(144);
pcolor( pfs{1}.adata.bins{:}, reshape(pfsa(:,u,8,6),pfs{1}.adata.binSizes') );
caxis   ([0,mrate]);
colormap('jet');
shading ('flat');
axis    ('xy');
xlim([-300,300]);
ylim([-300,300]);
title(num2str(units(u)));
waitforbuttonpress();
end
%%%>>>



%%%<<< PHASE dependent spatial shift information

% TRANSFORM Local Field Potential -> theta phase
Trial.lfp.filename = [Trial.name,'.lfp'];
phz = load(Trial,'lfp',thetaPhzChan).phase([6,12]);
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi + phzCorrection; % mv phzCorrection -> Trial prop
phz.data(phz.data<0) = phz.data(phz.data<0) + 2*pi;
phz.data(phz.data>2*pi) = phz.data(phz.data>2*pi) - 2*pi;
shifts = -100:20:100;

pfs = {};

for x = 1:numel(shifts),
    for y = 1:numel(shifts),
        for phase = 1:numel(binPhzc),
            
% CHECK existence of pfs object
            pargs.tag = ...
                ['egofield_refine_SW',num2str(pargs.SmoothingWeights(1)),                 ...
                 '_THP',num2str(phase),'-',num2str(numel(binPhzc)),                       ...
                 '_X',num2str(shifts(x)),'_Y',num2str(shifts(y))                          ...
                ];
            filepath = fullfile(Trial.spath, [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
        
            if exist(filepath,'file'),
                pfs{x,y,phase} = load(filepath,'Pfs').('Pfs');
                if overwrite,
                    pfs.purge_savefile();
                end;% if overwrite
            end;% if exist

        
            for unit = 1:numel(units),
                if unit==1 | electrode ~= spk.map(spk.map(:,1)==units(unit),2), % update phase state
                    pargs.spk = copy( spk );
                    electrode = 1;
                    pargs.states = copy( thetaState );
                    pargs.states.label = ['thetaPhz_',num2str(phase)];
                    pargs.states.data((phz(:,electrode) < binPhzs(phase) )                ...
                                      | (phz(:,electrode) >= binPhzs(phase+1)) ) = 0;
                    cast( pargs.states, 'TimePeriods' );
                    resInd = WithinRanges( pargs.spk.res, pargs.states.data );
                    pargs.spk.res = pargs.spk.res( resInd );
                    pargs.spk.clu = pargs.spk.clu( resInd );
                end;%if unit==1
            
                [mxr,mxp] = pft.maxRate(units(unit));
                pfsCenterHR = MTADfet.encapsulate(                                        ...
                    Trial,                                                                ... Trial
                    [multiprod(bsxfun(@minus,                                             ... data
                                      mxp+shifts([x,y]),                                  ...
                                      sq(xyz(:,'hcom',[1,2]))),                           ...
                               hvec,2,[2,3])],                                            ... 
                    sampleRate,                                                           ... SR
                    'egocentric_placefield',                                              ... name
                    'egopfs',                                                             ... label
                    'p'                                                                   ... key
                    );
                
                pargs.xyzp = pfsCenterHR;
                pargs.units  = units(unit);
                pfsArgs = struct2varargin(pargs);
                pfTemp = MTAApfs(pfTemp,pfsArgs{:});
                if unit==1,
                    try,
                        pfTemp.purge_savefile();
                    end;
                    pfTemp.save();        
                end;%if unit==1
            end;%for unit
            pfTemp.save();
            pfs{x,y} = pfTemp;    
            pfTemp = Trial;
        end;%for phase
    end;%for y
end;%for x
%%%>>>

pfsa = decapsulate_and_concatenate_mtaapfs(pfs(:)',repmat({units},[1,prod(size(pfs))]));
pfsa = reshape(pfsa,[size(pfsa,1),numel(units),size(pfs,1),size(pfs,2),size(pfs,3)]);

xyi = [21,1;21,11;21,21;11,1;11,11;11,21;1,1;1,11;1,21];
figure();
for i= 1:9,
subplot(3,3,i);
pcolor( pfs{1}.adata.bins{:}, reshape(pfsa(:,3,xyi(i,1),xyi(i,2)),pfs{1}.adata.binSizes') );
caxis   ([0,8.2]);
colormap('jet');
shading ('flat');
axis    ('xy');
end

bins = pfs{1}.adata.bins;
binSizes = pfs{1}.adata.binSizes;
width = binSizes(1);
height =binSizes(2);
radius = round(binSizes(1)/2)-find(bins{1}<=-300,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mazeMask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);

pfsavb = pfsa(logical(mazeMask(:)),:,:,:,:);

psi = permute(sum(repmat(1./sum(~isnan(pfsavb(:,:,1))),[size(pfsavb,1),1]) ...
              .*bsxfun(@rdivide,pfsavb,mean(pfsavb,'omitnan')) ...
              .*log2(bsxfun(@rdivide,pfsavb,mean(pfsavb,'omitnan'))),'omitnan'),...
              [3,4,5,2,1]);

figure();
nx = 4
ny = 5
for u = 1:21;
subplot2(ny,nx,1,1);        
plot(pft,units(u));
mpsi = LocalMinima2(-psi(:,:,2,u),0,3);
for p = 1:5,
subplot2(ny,nx,6-p,2);
    pcolor( pfs{1}.adata.bins{:}, reshape(pfsa(:,u,6,6,p),pfs{1}.adata.binSizes') );
    caxis   ([0,10]);
    colormap('jet');
    shading ('flat');
    axis    ('xy');
subplot2(ny,nx,6-p,3);
    hold('on');
    imagesc(shifts,shifts,psi(:,:,p,u)');
    plot(mpsi(2),mpsi(1),'*m');
    axis('xy');
    colorbar();
    caxis([0,1.8]);
subplot2(ny,nx,6-p,4);
    pcolor( pfs{1}.adata.bins{:}, reshape(pfsa(:,u,mpsi(1),mpsi(2),p),pfs{1}.adata.binSizes') );
    caxis   ([0,10]);
    colormap('jet');
    shading ('flat');
    axis    ('xy');
end;
title(num2str([u,units(u)]));
waitforbuttonpress();
end


figure();
%for u = 1:21;
    x = 6;
    y = 3;
subplot2(ny,nx,1,1);        
plot(pft,units(u));
for p = 1:5,
subplot2(ny,nx,6-p,2);
    pcolor( pfs{1}.adata.bins{:}, reshape(pfsa(:,u,6,6,p),pfs{1}.adata.binSizes') );
    caxis   ([0,15]);
    colormap('jet');
    shading ('flat');
    axis    ('xy');
subplot2(ny,nx,6-p,3);
    imagesc(psi(:,:,p,u));
    axis('xy');
    colorbar();
subplot2(ny,nx,6-p,4);
    pcolor( pfs{1}.adata.bins{:}, reshape(pfsa(:,u,x,y,p),pfs{1}.adata.binSizes') );
    caxis   ([0,15]);
    colormap('jet');
    shading ('flat');
    axis    ('xy');
end;
title(num2str(units(u)));
%waitforbuttonpress();
%end
