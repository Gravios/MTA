function Data = create(Data,Session,varargin)
%create(Data,Session,RefObj,spk,units,spikeWindow,mode,overwrite)
%
% Calculate the instantaneous firing rate of individual units
%
[RefObj,spk,units,spikeWindow,mode,overwrite] = DefaultArgs(varargin,{Session.lfp,[],[],0.05,'gauss',false});

if isa(RefObj,'MTAApfs')
    %nv units
    if isempty(units)
        units = cat(1,RefObj.data.clu);
    end
    Data.sampleRate = Session.xyz.sampleRate;
    ind = repmat({zeros(Session.xyz.size(1),1)},1,Session.xyz.size(3));
    for n = 1:numel(RefObj.adata.bins)
        [~,ind{n}] = NearestNeighbour(RefObj.adata.bins{n}',Session.xyz(:,Session.trackingMarker,n)');
    end
    Data.data = zeros(Session.xyz.size(1),numel(units));
    c = 1;
    for u = units
        rateMap = RefObj.data.rateMap(:,ismember(RefObj.data.clu,u),1);
        Data.data(:,c) = rateMap(sub2ind(RefObj.adata.binSizes',ind{:}));
        c = c+1;
    end
    
else
    
    Data.sampleRate = RefObj.sampleRate;

    % Load unit activity
    if isempty(spk)
        spk = Session.spk.copy();
        spk.create(Session,RefObj.sampleRate,[],units);
    end
    if spk.sampleRate ~= Data.sampleRate
        spk = spk.copy();
        spk.res = round(spk.res./spk.sampleRate.*Data.sampleRate);
        spk.sampleRate = Data.sampleRate;
    end
    
    
    if ~isempty(RefObj),
        dsize = RefObj.size(1);
    else
        dsize = diff(round(RefObj.sync.sync([1,end])*RefObj.sampleRate+1));
    end
    
    % Pre-allocate output array
    Data.data = zeros([dsize,numel(units)]);
    
    % If no units were provided load all units
    if isempty(units)
        units = 1:spk.map(end,1);
    end
    
    % create convolution window     
    swin = round(spikeWindow*RefObj.sampleRate);    

    %gwin = ones([swin,1]);
    switch mode
      case 'gauss'
        winLength = round(3*RefObj.sampleRate);
        gwin = gausswin(winLength,winLength/(spikeWindow*RefObj.sampleRate));
      case 'boxcar'
        gwin = ones([swin,1]);
        gwin = gwin./sum(gwin);
      case 'count'
        gwin = 1;
        spikeWindow = 1;
    end

    % accumulate and convolve activity of each unit
    for unit = units(:)'
        res = spk.res(spk.clu==unit);
        res(res==0) = 1;
        res(res>dsize) = dsize;
        Data.data(:,unit==units) = conv(accumarray(res,1,[dsize,1]),gwin,'same')+eps;
    end

    Data.spikeWindow = spikeWindow;
    Data.data(~nniz(RefObj),:) = 0;
    Data.origin = RefObj.origin;
    Data.sync = RefObj.sync;
end

