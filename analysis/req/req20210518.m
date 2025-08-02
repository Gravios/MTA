;; This buffer is for notes you don't want to save, and for Lisp evaluation.
;; If you want to create a file, visit that file with C-x C-f,
;; then enter the text in that file's own buffer.



MjgER2016_load_data();


cf(@(t,u) set_unit_set(t.spk,t,'pyramidal',u), Trials,units);
cf(@(t,u) set_unit_set(t.spk,t,'interneurons',u), Trials,unitsInts);

pyrs = cf(@(t) get_unit_set(t.spk,t,'pyramidal').units,    Trials);
ints = cf(@(t) get_unit_set(t.spk,t,'interneurons').units, Trials);

sampleRate = 250;


t = 20;
Trial = Trials{t};
units = pyrs{t};
xyz   = preproc_xyz(Trial,'trb');
xyz.resample(sampleRate);
rot = headRotation{t};

spk = Trial.load('spk',sampleRate,'walk+turn&theta',units);
pft = pfs_2d_theta(Trial,units);


Trial.lfp.filename = [Trial.name,'.lfp'];    
lfp = Trials{t}.load('lfp',sessionList(t).thetaRefGeneral);
phz = lfp.phase([5,13]);    
phz.data = unwrap(phz.data);
phz.resample(xyz);
phz.data = mod(phz.data+pi+phzCorrection(t),2*pi)-pi ; 
%phz.data(phz.data>pi) = phz.data(phz.data>pi)-2*pi;
%clear('lfp');



% COMPUTE head basis
hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
hvec = multiprod(hvec,                                  ...
                 [cos(rot),-sin(rot);sin(rot),cos(rot)],...
                 [2,3],                                 ...
                 [1,2]);



hvang = filter(copy(xyz),'ButFilter',4,2,'low');
xycoor = cat(2,...
             hvang(:,'spine_upper',[1,2])-hvang(:,'bcom',[1,2]),...
             hvang(:,'nose',[1,2])-hvang(:,'hcom',[1,2]));
hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
% Positive: CCW (Left)     Negative: CW (Right)
hvang.data = circ_dist(circshift(hvang.data(:,2),-10),...
                       circshift(hvang.data(:,2),+10));

hbang = filter(copy(xyz{t}),'ButFilter',3,30,'low');    
xycoor = cat(2,...
             hbang(:,'spine_upper',[1,2])-hbang(:,'bcom',[1,2]),...
             hbang(:,'nose',[1,2])-hbang(:,'hcom',[1,2]));
hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
hbang.data = circ_dist(hbang.data(:,2),hbang.data(:,1));
hbang.data = hbang.data - 0.25;


hrvfl = fet_href_HXY(Trial,sampleRate,false,'trb',4);



% GET theta state behaviors, minus rear
thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);





pfTemp = Trial;

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = units;
pargs.tag          = 'egofield';
pargs.binDims      = [20, 20]; % X Y 
pargs.SmoothingWeights = [3, 3];
pargs.halfsample   = false;
pargs.numIter      = 1;   
pargs.boundaryLimits = [-410,410;-410,410];
pargs.states       = '';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
electrode = 0;

% CHECK existence of pfs object
if pargs.SmoothingWeights(1)~=2,
    pargs.tag = ['egofield_SW',num2str(pargs.SmoothingWeights(1))];
else
    pargs.tag = 'egofield';
end
filepath = fullfile(Trial.spath,...
                    [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
    
if exist(filepath,'file'),
    pfs = load(filepath);
    pfs = pfs.Pfs;
    if overwrite & isa(pfs,'MTAApfs'),
        pfs.purge_savefile();
    else,
        return;
    end;% if overwrite
end;% if exist




figure()
hold('on');
plot(-hbang(res));
plot(-hvang(res));
plot(hrvfl(res,2)/50);
plot(,'k')



hfig = figure();
% $$$ mxp = [0,0];

%mxp = [126,-180];
clf();
unit = 21;
[~,mxp] = pft.maxRate(unit);
c = 100;
res = spk(unit);
res = res(hrvfl(res,1)>5);
sphz = phz(res);
sxyz = sq(xyz(res,'hcom',[1,2]));
shbang = -hbang(res);
%sbang = (cos(sphz+pi/2)+1)/2.*(-hbang(res));
sbang = (cos(sphz+pi/2)+1)/2.*(-hbang(res)+hrvfl(res,2)/40);
%sbang = (-hbang(res)+hrvfl(res,2)/60);
%sbang = sign(sphz).*(-hbang(res)-hvang(res)+hrvfl(res,2)/40);
%ascale = -1;
%ascale = -0.75;
ascale = -0.25;
%ascale = 0.75;
ascale = -0.5;
ascale = 0.001;

srot = cat(3,[cos(sbang.*ascale),-sin(sbang.*ascale)],...
             [sin(sbang.*ascale), cos(sbang.*ascale)]);
shvec = hvec(res,:,:);
spkEgoXY = multiprod(bsxfun(@minus,                  ...
                            mxp,                     ...
                            sxyz),...
                     shvec,2,[2,3]);
sax  = gobjects([1,2]);
pax  = gobjects([1,2]);
sax(1) = subplot2(3,2,1,1);
    hold(sax(1),'on');
    plot(pft,unit,1,'colorbar');
    pax(1) = scatter(mxp(1),mxp(2),20,'r','Filled');
sax(2) = subplot2(3,2,2,1);
    hold(sax(2),'on');
    pind = sphz>0;    
    %plot(xyz(thetaState,'hcom',1),xyz(thetaState,'hcom',2)
    pax(2) = scatter(spkEgoXY(pind,1),spkEgoXY(pind,2),10,sphz(pind),'Filled');
    colormap(gca(),'hsv');
    colorbar();
    caxis([-pi,pi]);    
    xlim(sax(2),[-300,300])    
    ylim(sax(2),[-300,300])        
% $$$ while true
    if ~ascale
        waitforbuttonpress();
        mxp = get(hfig.CurrentAxes,'CurrentPoint');
        mxp = mxp(1,1:2);
        pax(1).XData = mxp(1);
        pax(1).YData = mxp(2);
    end    
    spkEgoXY_down = multiprod(...
        multiprod(bsxfun(@minus,                  ...
                         mxp,                     ...
                         sxyz(pind,:)),...
                  shvec(pind,:,:),2,[2,3]),...
        srot(pind,:,:),...
        [2],...
        [2,3]);
    pax(2).XData = spkEgoXY_down(:,1);
    pax(2).YData = spkEgoXY_down(:,2);
subplot2(3,2,2,2);
    hist2([spkEgoXY_down(:,1),spkEgoXY_down(:,2)],...
          linspace(-300,300,11),...
          linspace(-300,300,11));
    colormap(gca(),'jet');
    colorbar();
    caxis([0,c]);
    
sax(3) = subplot2(3,2,3,1);
    hold(sax(3),'on');
    pind = sphz<0;
    %plot(xyz(thetaState,'hcom',1),xyz(thetaState,'hcom',2)
    pax(3) = scatter(spkEgoXY(pind,1),spkEgoXY(pind,2),10,sphz(pind),'Filled');
    %pax(3) = scatter(spkEgoXY(pind,1),spkEgoXY(pind,2),10,shbang(pind),'Filled');    
    caxis([-pi,pi]);
    colormap(gca(),'hsv');
    colorbar();
    xlim(sax(3),[-300,300])    
    ylim(sax(3),[-300,300])        
    spkEgoXY_up = multiprod(...
        multiprod(bsxfun(@minus,                  ...
                         mxp,                     ...
                         sxyz(pind,:)),...
                  shvec(pind,:,:),2,[2,3]),...
        srot(pind,:,:),...
        [2],...
        [2,3]);
    pax(3).XData = spkEgoXY_up(:,1);
    pax(3).YData = spkEgoXY_up(:,2);
subplot2(3,2,3,2);
    hist2([spkEgoXY_up(:,1),spkEgoXY_up(:,2)],...
          linspace(-300,300,11),...
          linspace(-300,300,11));
    colormap(gca(),'jet');
    colorbar();
    caxis([0,c]);
    
    %end













unit = 21;
[~,mxp] = pft.maxRate(unit);
res = spk(unit);
sphz = phz(res);
sxyz = sq(xyz(res,'hcom',[1,2]));
shbang = -hbang(res);

shvec = hvec(res,:,:);
sbang = (cos(sphz+pi/2)+1)/2.*(-hbang(res)+hrvfl(res,2)/60);
sbang = sbang(randperm(numel(sbang)));

msize = [21,21];
ndims = numel(msize);
sind = cell(1,ndims);
for i = 1:ndims
    sind{i} = linspace(-round(msize(i)/2),round(msize(i)/2),msize(i));
end
[sind{:}] = ndgrid(sind{:});


SmoothingWeights = [1.1,1.1];
for i = 1:ndims,
    sind{i} = sind{i}.^2/SmoothingWeights(i)^2/2;
end
Smoother = exp(sum(-cat(ndims+1,sind{:}),ndims+1));
Smoother = Smoother./sum(Smoother(:));

% SMOOTH occupancy
% SMOOTH Spike count

si = [];
for ascale = -2:0.1:2
srot = cat(3,[cos(sbang.*ascale),-sin(sbang.*ascale)],...
             [sin(sbang.*ascale), cos(sbang.*ascale)]);
pind = sphz<0;
spkEgoXY_up = multiprod(...
    multiprod(bsxfun(@minus,                  ...
                     mxp,                     ...
                     sxyz(pind,:)),...
              shvec(pind,:,:),2,[2,3]),...
    srot(pind,:,:),...
    [2],...
    [2,3]);
out =hist2([spkEgoXY_up(:,1),spkEgoXY_up(:,2)],...
           linspace(-300,300,22),...
           linspace(-300,300,22));
    out   = convn(out, Smoother,'same');

si(end+1) = sum(1./numel(out).*out(:)./mean(out(:)).*log2(out(:)./mean(out(:))),'omitnan');
end
figure();
hold('on');
plot(-2:0.1:2,si)







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
    end;% if
    
    [mxr,mxp] = pft.maxRate(units(unit));
    pfsCenterHR = MTADfet.encapsulate(Trial,                                            ...
                                      [bsxfun(@plus,                                    ...
                                              multiprod(bsxfun(@minus,                  ...
                                                               mxp,                     ...
                                                               sq(xyz(:,'hcom',[1,2]))),...
                                                         hvec,2,[2,3]),                 ...
                                              headCenterCorrection)],                   ...
                                       sampleRate,                                      ...
                                      'egocentric_placefield',                          ...
                                      'egopfs',                                         ...
                                      'p'                                               ...
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
    end;% if 
end;% for unit
pfTemp.save();
pfs = pfTemp;    
pfTemp = Trial;

