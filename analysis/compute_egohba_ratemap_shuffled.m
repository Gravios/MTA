function [pfs] = compute_egohba_ratemap_shuffled(Trial,units,xyz,spk,pft,overwrite)


sampleRate = xyz.sampleRate;
binPhzs = linspace(0.5,2*pi-0.5,4);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;
binHbas = [-1.2,-0.2,0.2,1.2];
binHbac = (binHbas(1:end-1)+binHbas(2:end))./2;
pfs = cell([numel(binPhzc)]);
verbose = true;

if verbose,
    disp(['[status]        compute_egohba_ratemap: processing trial: ',Trial.filebase]);
end


if isempty(units),
    return;
end;% if


if overwrite,
    hba = fet_head_body_angle(Trial, 'xyz', xyz);
    phz = load_theta_phase(Trial, xyz.sampleRate);
    rot = transform_vector_to_rotation_matrix( ...
             xyz,{'hcom','nose'}, Trial.meta.correction.headYaw);

    % GET theta state behaviors, minus rear
    thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);
end

pfTemp = Trial;

pargs = get_default_args('MjgER2016','MTAApfs','struct');        
pargs.units        = units;
pargs.tag          = 'egofield';
pargs.binDims      = [20, 20, 0.6];                           % X Y HBA
pargs.SmoothingWeights = [3, 3, 0.4];                     % X Y HBA
pargs.type         = 'xyw';
pargs.spkShuffle   = false;
pargs.posShuffle   = true;
pargs.halfsample   = false;
pargs.numIter      = 100;   
pargs.boundaryLimits = [-410,410;-410,410;-1.5,1.5];
pargs.states       = '';
pargs.overwrite    = true;
pargs.autoSaveFlag = false;    
pargs.posShuffleDims = 3;
electrode = 0;

% Don't judge me
if pargs.SmoothingWeights(1)~=2,
    stag = ['_SW',num2str(pargs.SmoothingWeights(1))];
else
    stag = '';
end

for phase = 1:numel(binPhzc)
% CHECK existence of pfs object
    pargs.tag = ['egofield_theta_phase_hba_shuffled_',num2str(phase),stag];
    filepath = fullfile(Trial.spath, [Trial.filebase,'.Pfs.',pargs.tag,'.mat']);

    
    if exist(filepath,'file'),
        pfs{phase} = load(filepath).Pfs;
        if overwrite,
            pfs{phase}.purge_savefile();
        else,
            continue;
        end;% if
    end;% if exist
    
        
    for unit = 1:numel(units),
        if unit==1 | electrode~=spk.map(spk.map(:,1)==units(unit),2), % update phase state            
            pargs.spk = copy(spk);
            electrode = 1;
            %electrode = spk.map(spk.map(:,1)==units(unit),2);
            pargs.states = copy(thetaState);
            pargs.states.label = ['thetaPhz_',num2str(phase)];
            pargs.states.data((phz(:,electrode) < binPhzs(phase) )    ...
                              | (phz(:,electrode) >= binPhzs(phase+1)) ) = 0;
            cast(pargs.states,'TimePeriods');
            resInd = WithinRanges(pargs.spk.res,pargs.states.data);
            pargs.spk.res = pargs.spk.res(resInd);
            pargs.spk.clu = pargs.spk.clu(resInd);
        end;% if
        
        [mxr,mxp] = pft.maxRate(units(unit));
        pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
                                          [bsxfun(                                         ...
                                              @plus,                                       ...
                                              multiprod(bsxfun(@minus,                     ...
                                                          mxp,                             ...
                                                          sq(xyz(:,'hcom',[1,2]))),        ...
                                                     rot,2,[2,3]),                         ...
                                              Trial.meta.correction.headCenter),           ...
                                           hba.data],                                      ...
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
            try
                pfTemp.purge_savefile();
            end
            pfTemp.save();        
        end% if 
    end;% for unit
    pfTemp.save();
    pfs{phase} = pfTemp;    
    pfTemp = Trial;
end;% for phase




%% PLOTING EXAMPLES
% $$$ lims = {[-250,250],[-250,250]};
% $$$ figure();
% $$$ al = 1:5;
% $$$ numAng = numel(al);
% $$$ sax = gobjects([0,1]);
% $$$ for p = 4
% $$$     rmap = plot(pfs{p},15);
% $$$     for a = 1:numAng,
% $$$         sax(end+1) = subplot2(1,5,1,a);
% $$$         pcolor(pfs{p}.adata.bins{1},...
% $$$                        pfs{p}.adata.bins{2},...
% $$$                        rmap(:,:,al(a)));
% $$$         
% $$$         caxis   (sax(end),[0,10]);
% $$$         colormap(sax(end),'jet');
% $$$         shading (sax(end),'flat');
% $$$         axis    (sax(end),'xy');
% $$$         xlim    (sax(end),lims{1});
% $$$         ylim    (sax(end),lims{2});        
% $$$         
% $$$         Lines([],0,'k');
% $$$         Lines(0,[],'k');
% $$$         
% $$$         set(sax(end),'XTick',[]);
% $$$         set(sax(end),'YTick',[]);        
% $$$         
% $$$         % ADD subject
% $$$ % $$$         if p %== 4,
% $$$ % $$$             subject = struct(rat);
% $$$ % $$$             subject = update_subject_patch(subject,'head',[], false,hbaBinEdg,hbaBinCtr);
% $$$ % $$$             subject = update_subject_patch(subject,'body', numAng+1-a,  true,hbaBinEdg,hbaBinCtr);
% $$$ % $$$             patch(subject.body.patch.vert{:},   [0.75,0.75,0.75]);
% $$$ % $$$             patch(subject.head.patch.vert{:},   [0.75,0.75,0.75]);
% $$$ % $$$             patch(subject.body.overlay.vert{:},[0.75,0.50,0.50],'FaceAlpha',0.3);
% $$$ % $$$             line(subject.head.midline.vert{:}, 'Color', subject.head.midline.color);
% $$$ % $$$             line(subject.body.midline.vert{:}, 'Color', subject.body.midline.color);
% $$$ % $$$         end
% $$$     end
% $$$ end
