function Data = create(Data,Session,varargin)
%create(Data,Session,DataObj,state,units,twin,overwrite)
%
% Calculate the instantaneous firing rate of individual units
%
[DataObj,state,units,twin,overwrite,mode] = DefaultArgs(varargin,{Session.lfp,[],[],0.05,false,'gauss'});

if isa(DataObj,'MTAApfs'),
    %nv units
    if isempty(units),
        units = cat(1,DataObj.data.clu);
    end
    Data.sampleRate = Session.xyz.sampleRate;
    ind = repmat({zeros(Session.xyz.size(1),1)},1,Session.xyz.size(3));
    for n = 1:numel(DataObj.adata.bins),
        [~,ind{n}] = NearestNeighbour(DataObj.adata.bins{n}',Session.xyz(:,Session.trackingMarker,n)');
    end
    Data.data = zeros(Session.xyz.size(1),numel(units));
    c = 1;
    for u = units,
        rateMap = DataObj.data.rateMap(:,ismember(DataObj.data.clu,u),1);
        Data.data(:,c) = rateMap(sub2ind(DataObj.adata.binSizes',ind{:}));
        c = c+1;
    end
    
else
    
    Data.sampleRate = DataObj.sampleRate;

    % Load unit activity
    spk = Session.spk.copy();
    spk.create(Session,DataObj.sampleRate,state,units);

    
    if ~isempty(DataObj),
        dsize = DataObj.size(1);
    else
        dsize = diff(round(DataObj.sync.sync([1,end])*DataObj.sampleRate+1));
    end
    
    % Pre-allocate output array
    Data.data = zeros([dsize,numel(units)]);
    
    % If no units were provided load all units
    if isempty(units), units = 1:spk.map(end,1);end
    
    % create convolution window     
    swin = round(twin*DataObj.sampleRate);    

    %gwin = ones([swin,1]);
    switch mode
        case 'gauss'
        winLength = round(3*DataObj.sampleRate);
        gwin = gausswin(winLength,winLength/(twin*DataObj.sampleRate));
      case 'boxcar'
        gwin = ones([swin,1]);
        gwin = gwin./sum(gwin);
      case 'count'
        gwin = 1;
        twin = 1;
    end

    % accumulate and convolve activity of each unit
    for unit = units(:)'
        res = spk.res(spk.clu==unit);
        res(res==0) = 1;
        res(res>dsize) = dsize;
        Data.data(:,unit==units) = conv(accumarray(res,1,[dsize,1]),gwin,'same')+eps;
    end

    Data.data(~nniz(DataObj),:) = 0;
    Data.origin = DataObj.origin;
    Data.sync = DataObj.sync;
end

