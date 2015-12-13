function phs = phase(Data,varargin)
%function phs = phase(Data,varargin)
%[freq_range,n] = DefaultArgs(varargin,{[6,12],3});
[freq_range,n] = DefaultArgs(varargin,{[6,12],3});
tbp = Data.copy;
tbp.filter('ButFilter',n,freq_range,'bandpass');
tbp_hilbert = Shilbert(Data.data);
tbp_phase = phase(tbp_hilbert);
DataClass = class(Data);
phs = feval(DataClass,'data',tbp_phase,...
            'sampleRate',Data.sampleRate,...
            'sync',Data.sync.copy,...
            'origin',Data.origin);
end
