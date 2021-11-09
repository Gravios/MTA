function phs = phase(Data,varargin)
%function phs = phase(Data,varargin)
%[freq_range,n] = DefaultArgs(varargin,{[6,12],3});
[freq_range,n] = DefaultArgs(varargin,{[6,12],3});

tbp = Data.copy;
tbp.filter('ButFilter',n,freq_range,'bandpass');

tbpHilbert = zeros(size(Data));
tbpHilbert(nniz(Data),:) = hilbert(tbp.data(nniz(Data),:));

DataClass = class(Data);
phs = feval(DataClass,...
                  'data',phase(tbpHilbert),...
            'sampleRate',Data.sampleRate,...
                  'sync',Data.sync.copy,...
                'origin',Data.origin);
end
