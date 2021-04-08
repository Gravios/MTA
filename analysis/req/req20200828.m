%req20200828
% sanity check
% find best parameters for placefield estimation
%
% Method: kernel density estimation
% Similar to: overdispertion
%     Difference: use placefield estimate instead of inhomogeneous poisson process model
%
% compute the expected number of spikes based on the ratemap
%     integrate the expected number of spikes along the trajectory weighted by the sampling rate.
%     integrate the observed number of spikes along the trajectory
%
% integrate along MTADufr 



%configure_default_args();

Trial = MTATrial.validate('jg05-20120317.cof.all');

    
pf = {};

% loop: bin size 10:10:100
binSizes = 10:10:100;
% loop: smoothing weights
sigma = 0.5:0.25:3;
for b = 1:numel(binSizes),
    for s = 1:numel(sigma),
        pf{b,s} = compute_ratemaps(Trial,...
                                   'pfsArgs', struct('states',           'theta-groom-sit',            ...
                                                     'binDims',          binSizes(b).*[1,1],           ...
                                                     'SmoothingWeights', sigma(s).*[1,1],              ...
                                                     'numIter',          1,                            ...
                                                     'boundaryLimits',   [-500,500;-500,500],          ...
                                                     'halfsample',       false),...
                                   'overwrite',true);
        %---------------------------------------------------------------------------------------------------
    end
end


figure();
for b = 1:numel(binSizes),
    for s = 1:numel(sigma),
        subplot2(numel(binSizes),numel(sigma),b,s);
        plot(pf{b,s},10);
    end
end


xyz = preproc_xyz(Trial,'trb');
xyz.resample(16);


ufr = Trial.load('ufr',xyz,[],[],[],[],'count');

binLims = [-500,500;-500,500]




% select 500ms trajectory segments 
%     where -
%        theta state


tx = sq(xyz(Trial.stc{'t'},'hcom',[1,2]));
tu = sq(ufr(Trial.stc{'t'},:));

nind = nniz(tx);

% get the r(x,y) from ratemaps



% CONSTRUCT map of expected rate given position
win = 16;
sse = zeros(numel(binSizes),numel(sigma),size(tu,2));

featureCA = mat2cell(tx(nind,:),sum(nind),ones([1,size(tx,2)]));

for b = 1:numel(binSizes),
    for s = 1:numel(sigma),
        tr = zeros(size(tx,1),size(tu,2));                    
        for u = 1:numel(units),
            rateMap = pf{b,s}.plot(u,'mean',[],[],false,[],false,[]);
            rateMap(isnan(rateMap)) = 0;
            tr(nind,u) = interpn(pf{b,s}.adata.bins{:},rateMap,featureCA{:},'linear');    
        end
        tr(tr<0) = 0;
        trs = sq(sum(GetSegs(tr,1:win:size(tr,1),win,nan),'omitnan')).*/win;
        tus = sq(sum(GetSegs(tu,1:win:size(tu,1),win,nan),'omitnan'));
        sse(b,s,:) = sum((trs-tus).^2,'omitnan')./size(trs,1);
    end
end

figure()
    imagesc(sse(:,:,10)');

cmap = cool(11);
figure();
    hold('on');
    for c = 1:10,
        plot(sse(c,:,10),'Color',cmap(c,:));
    end
    
figure();
hold('on');
plot(diff(sse(2,:,10)))



figure();
for u = 1:size(Trial.spk.map,1),
    plot(pf{2,5},u,1,'colorbar');
    title(u);
    waitforbuttonpress();
    clf();
end

units = select_placefields(Trial);


% good:  7,10,21
% bad : 20,26
% ugly:  6, 7,12,19,20,22,41
% 2low:  8,21
% 

states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
          'pause&theta','lpause&theta','hpause&theta'};
pfs = pfs_2d_states(Trial,1:76,[],states,false,[],true,struct('numIter',1,'halfsample',false));


figure();
for u = 51:size(Trial.spk.map,1),
    mxr = max(cell2mat(cf(@(p) p.maxRate(u),pfs)));
    subplot(1,8,1)
    plot(pf{2,5},u,1,'colorbar',[0,mxr],'colorMap',@jet);
    title(u);   
    for s = 1:7,
        subplot(1,8,s+1);
        plot(pfs{s},u,1,'colorbar',[0,mxr],'colorMap',@jet);
        title(states{s});
    end
    waitforbuttonpress();
    clf();
end
