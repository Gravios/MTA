
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
thetaChan = sessionList(20).thetaRefGeneral;

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

for x = 1:11%numel(shifts),
    for y = 1:numel(shifts),
        pargs.tag = ['egofield_refine_SW',num2str(pargs.SmoothingWeights(1)),...
                     '_X',num2str(shifts(x)),'_Y',num2str(shifts(y))];
        filepath = fullfile(Trial.spath, [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);
            
        if exist(filepath,'file'),
            pfs{x,y} = load(filepath);
            if overwrite,
                pfs.purge_savefile();
            end;% if overwrite
        end;% if exist
    end
end

        
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



pfsa = decapsulate_and_concatenate_mtaapfs(pfs(:)',repmat({units},[1,prod(size(pfs))]));

bsi = sum(repmat(1./sum(~isnan(brmaps(:,:,1))),[size(brmaps,1),1]).*bsxfun(@rdivide,brmaps,mean(brmaps,'omitnan')) ...
          .*log2(bsxfun(@rdivide,rmaps,mean(brmaps,'omitnan'))),'omitnan');
