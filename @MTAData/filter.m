function Data = filter(Data,mode,varargin)
% Data = filter(Data,varargin)
% Data.data = reshape(Filter0(win,Data.data),Data.size);
%
% TODO Adding ButFilter
switch Data.type
  case 'TimeSeries'
    switch mode
      case 'gauss'
        [Window] = DefaultArgs(varargin,{.05});
        if numel(Window)==1,
            Window = round(Window.*Data.sampleRate);
            Window = Window+double(mod(Window,2)==0);
            Window = gausswin(Window)./sum(gausswin(Window));
        end
        zind = nniz(Data);
        Data.data = reshape(Filter0(Window,Data.data),Data.size);
        Data.data(zind) = 0;
      case 'ButFilter'
        [order,freq,flag] = DefaultArgs(varargin,{3,4,'low'});
        freq = freq/(Data.sampleRate/2);
        nind = ~nniz(Data);
        nd = [];
        if sum(nind)>0,
            nd = Data.data(nind,:,:,:,:);
            Data.data(nind,:,:,:,:) = 0;
        end
        Data.data = ButFilter(Data.data,order,freq,flag);
        if ~isempty(nd),
            Data.data(nind,:,:,:,:) = nd;
        end
      case 'RectFilter'
        Data.data = RectFilter(Data.data,varargin{:});
        
      otherwise
        error('MTA:MTAData:filter:ModeNotFound');
    end
  otherwise
    return
end
end

