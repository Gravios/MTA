function [fet,varargout] = fet_rds(Trial,varargin)
% [rhm,fs,ts] = fet_rhm(Trial,varargin)
% [sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'spectral',1});
% Need to update the spectral window size to adapt to the xyz.sampleRate
%

Trial = MTATrial.validate(Trial);

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate',             Trial.xyz.sampleRate,                              ...
                 'normalize'    ,             false,                                             ...
                 'procOpts'     ,             {{'SPLINE_SPINE_HEAD_EQD'}},                       ...
                 'referenceTrial',            'Ed05-20140529.ont.all',                           ...
                 'referenceFeature',          ''                                                 ...
);
[newSampleRate,normalize,procOpts,referenceTrial,referenceFeature] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

xyz = preproc_xyz(Trial,procOpts);
xyz.filter('RectFilter');
pch = fet_HB_pitchB(Trial,varargin);


mxy = decompose_xy_motion_wrt_body(Trial,'maxNumComponents',2);
mxy.resample(xyz);

% TAG creation -------------------------------------------------------------------------------------
hash = DataHash(struct('sampleRate',                   newSampleRate,                            ...
                       'normalize',                    normalize,                                ...
                       'xyz',                          xyz.hash,                                 ...
                       'pch',                          pch.hash,                                 ...
                       'mxy',                          pch.hash                                  ...
));
%---------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'reduced feature set','fet_rds','r');

fet.data = nunity([pch.data,...
                   mxy(:,1),...
                   diff([circshift(xyz(:,'hcom',3),-1),circshift(xyz(:,'hcom',3),1)],1,2)]);

% END MAIN -----------------------------------------------------------------------------------------