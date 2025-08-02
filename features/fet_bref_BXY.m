function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_bref_BXY(Trial,varargin)
% function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_bref_BXY(Trial,varargin)
% 
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%     filterCutoff:  numeric,  (4)                    - low pass filter boundary frequency
%     theta:         numeric,  (0)                    - rotation of the basis in radians
%
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate'  , Trial.xyz.sampleRate,                                        ...
                 'normalize'      , false,                                                       ...
                 'procOpts'       , {{''}},                                                      ...
                 'filterCutoff'   , 4,                                                           ...
                 'theta'          , 0                                                            ...
);
[newSampleRate,normalize,procOpts,filterCutoff,theta] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'head referenced position and motion','fet_href_H','h');

% PREPROC xyz
% FILTER xyz
xyz = preproc_xyz(Trial,procOpts);
xyz.resample(newSampleRate);
xyz.data(~nniz(xyz),:,:) = 0;
xyz.filter('ButFilter',4,filterCutoff,'low');
bcom = xyz(:,'bcom',[1,2]);

% GENERATE orthogonal basis, origin: head's center of mass
nvec = xyz(:,'spine_upper',[1,2])-bcom;
nvec = bsxfun(@rdivide,nvec,sqrt(sum((nvec).^2,3))); 
nvec = cat(2,nvec,permute(sq(nvec)*[0,-1;1,0],[1,3,2]));
nvec =  multiprod([cos(theta),-sin(theta);sin(theta),cos(theta)],nvec,[1,2],[2,3]);

% COMPUTE hcom projection onto head reference
uvec = circshift(bcom,-1)-circshift(bcom,1);
fet.data(:,1) = dot(uvec,nvec(:,1,:),3).*(xyz.sampleRate/10)*0.5;
fet.data(:,2) = dot(uvec,nvec(:,2,:),3).*(xyz.sampleRate/10)*0.5;

% DIAGNOSTIC figure
% $$$ figure,
% $$$ ind = Trial.stc{'w'};
% $$$ hist2(fet(ind,:),...
% $$$       linspace(-100,100,100),...
% $$$       linspace(-100,100,100));

featureTitles = {};
featureDesc = {};
if nargout>1,
end


%---------------------------------------------------------------------------------------------------





