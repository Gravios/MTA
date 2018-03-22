%
o
Trial = MTATrial.validate('jg05-20120317.cof.all');
Trial.load('stc','msnn_ppsvd_raux');
xyz = preproc_xyz(Trial,'trb');
units = select_placefields(Trial);


pargs.units = units;
pargs.states = 'loc&theta';
pargs.overwrite = true;
pargs.numIter = 101;

pfs = {};
sampleRates = [2.^[1:7]];
for r = 1:numel(sampleRates),
    pargs.xyzp = xyz.copy('empty');
    pargs.xyzp.data = sq(xyz(:,'nose',[1,2]));
    pargs.xyzp.resample(sampleRates(r));
    pfsArgs = struct2varargin(pargs);
    pfs{r} = MTAApfs(Trial,pfsArgs{:});
end


figure();
for u = 1:numel(units),    
    clf();
    for r = 1:numel(sampleRates),
        subplot(1,numel(sampleRates),r);
        plot(pfs{r},units(u),'mean',true,[],false,0.99);
        title(num2str(sampleRates(r)));
    end
    waitforbuttonpress();
end



pargs.numIter = 10001;
sampleRates = 16;
pargs.xyzp = xyz.copy('empty');
pargs.xyzp.data = sq(xyz(:,'nose',[1,2]));
pargs.xyzp.resample(16);
pfsArgs = struct2varargin(pargs);
pff = MTAApfs(Trial,pfsArgs{:});

pargs.numIter = 10001;
pargs.states = 'theta-groom-sit';
pargs.xyzp = xyz.copy('empty');
pargs.xyzp.data = sq(xyz(:,'nose',[1,2]));
pargs.xyzp.resample(16);
pfsArgs = struct2varargin(pargs);
pft = MTAApfs(Trial,pfsArgs{:});

pargs.numIter = 1001;
pargs.states = 'theta-groom-sit';
pargs.overwrite = true;
pargs.xyzp = xyz.copy('empty');
pargs.xyzp.data = sq(xyz(:,'nose',[1,2]));
pargs.xyzp.resample(16);
pfsArgs = struct2varargin(pargs);
pftn = MTAApfs(Trial,pfsArgs{:});



figure();

for u = 1:numel(units),
    clf();
    subplot(3,4,1);
    plot(pff,units(u),'mean',true,[],true,0.5);
    subplot(3,4,2);
    plot(pff,units(u),'std', true,[],true,0.5);
    subplot(3,4,3);
    plot(pff,units(u),'snr', true,[],true,0.5);
    subplot(3,4,4);
    plot(pff,units(u),'snrs', true,[],true,0.5);

    subplot(3,4,5);
    plot(pft,units(u),'mean',true,[],true,0.5);
    subplot(3,4,6);
    plot(pft,units(u),'std', true,[],true,0.5);
    subplot(3,4,7);
    plot(pft,units(u),'snr', true,[],true,0.5);
    subplot(3,4,8);
    plot(pft,units(u),'snrs', true,[],true,0.5);

    subplot(3,4,9);
    plot(pftn,units(u),'mean',true,[],true,0.5);
    subplot(3,4,10);
    plot(pftn,units(u),'std', true,[],true,0.5);
    subplot(3,4,11);
    plot(pftn,units(u),'snr', true,[],true,0.5);
    subplot(3,4,12);
    plot(pftn,units(u),'snrs', true,[],true,0.5);

    waitforbuttonpress();
end



pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.numIter          = 101;
pargs.binDims          = [50,50];
pargs.SmoothingWeights =  [0.75 0.75];
pargs.states = 'theta-groom-sit';
pargs.overwrite = true;
pargs.xyzp = xyz.copy('empty');
pargs.xyzp.data = sq(xyz(:,'nose',[1,2]));
pargs.xyzp.resample(16);
pfsArgs = struct2varargin(pargs);
pftg = MTAApfs(Trial,pfsArgs{:});


4*(50/200).^2

figure();
ny = 4;
hsThresh = 0.9;
iCirc = false;
iColorbar = true;
for u = 1:numel(units),
    clf();
    subplot(ny,4,1);
    plot(pftn,units(u),'mean',iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,2);
    plot(pftn,units(u),'std', iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,3);
    plot(pftn,units(u),'snr', iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,4);
    plot(pftn,units(u),'snrs',iColorbar,[],iCirc,hsThresh);
    
    l = 4;
    subplot(ny,4,l+1);
    plot(pft,units(u),'mean',iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+2);
    plot(pft,units(u),'std', iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+3);
    plot(pft,units(u),'snr', iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+4);
    plot(pft,units(u),'snrs',iColorbar,[],iCirc,hsThresh);

    l = 8;
    subplot(ny,4,l+1);
    plot(pftb,units(u),'mean',iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+2);
    plot(pftb,units(u),'std', iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+3);
    plot(pftb,units(u),'snr', iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+4);
    plot(pftb,units(u),'snrs',iColorbar,[],iCirc,hsThresh);

    l = 12;
    subplot(ny,4,l+1);
    plot(pftg,units(u),'mean',iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+2);
    plot(pftg,units(u),'std', iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+3);
    plot(pftg,units(u),'snr', iColorbar,[],iCirc,hsThresh);
    subplot(ny,4,l+4);
    plot(pftg,units(u),'snrs',iColorbar,[],iCirc,hsThresh);
    
    waitforbuttonpress();
end

sind = cell([1,2]);
[sind{:}] = ndgrid(linspace(-500,500,100),linspace(-500,500,100));

bins = linspace(-500,500,100);

width = 100;
height = 100;
radius = round(50)-find(bins<-450,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
pmask = double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius);
pmask(pmask==0)=nan;




rmap = plot(pftg,94,'mean',false,[],true);
rmask = double(isnan(rmap));
rmap(isnan(rmap)) = 0;

nmap = interp2(pftg.adata.bins{:},rmap,sind{:},'cubic');
nmask = interp2(pftg.adata.bins{:},rmask,sind{:},'linear');
rmap(rmask==1) = nan;
nmap(nmask>=0.75) = nan;
nmap(nmap<0) = 0;

figure,
subplot(231);
imagescnan(rot90(rmap))
subplot(232);
imagescnan(rot90(nmap'))
subplot(233);
plot(pftg,94,'mean',true,[],true);

rmap = plot(pftg,units(21),'mean',false,[],true,0.9);
rmask = double(isnan(rmap));
rmap(isnan(rmap)) = 0;
nmap = interp2(pftg.adata.bins{:},rmap,sind{:},'cubic');
nmask = interp2(pftg.adata.bins{:},rmask,sind{:},'linear');
rmap(rmask==1) = nan;
nmap(nmask>=0.75) = nan;
nmap(nmap<0) = 0;


subplot(234);
imagescnan(rot90(rmap))
subplot(235);
imagescnan(rot90(nmap').*pmask)
subplot(236);
plot(pftg,units(21),'mean',iColorbar,[],true,0.99);
