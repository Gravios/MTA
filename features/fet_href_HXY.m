function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_href_HXY(Trial,varargin)
% function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_href_HXY(Trial,varargin)
% 
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate'  , Trial.xyz.sampleRate,                                        ...
                 'normalize'      , false,                                                       ...
                 'procOpts'       , {{''}},                                                      ...
                 'filterCutoff'   , 4,                                                         ...
                 'theta'          , 0 ,                                                          ...
                 'markers'        , {{'head_back','head_left','head_front','head_right'}}        ...
);
[newSampleRate,normalize,procOpts,filterCutoff,theta,markers] =                               ...
    DefaultArgs(varargin,defargs,'--struct');
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

hcom = xyz(:,'hcom',[1,2]);

% GENERATE orthogonal basis, origin: head's center of mass
nvec = xyz(:,'nose',[1,2])-hcom;
nvec = bsxfun(@rdivide,nvec,sqrt(sum((nvec).^2,3))); 
nvec = cat(2,nvec,permute(sq(nvec)*[0,-1;1,0],[1,3,2]));
nvec =  multiprod([cos(theta),-sin(theta);sin(theta),cos(theta)],nvec,[1,2],[2,3]);

% COMPUTE hcom projection onto head reference
uvec = circshift(hcom,-1)-circshift(hcom,1);
fet.data(:,1) = dot(uvec,nvec(:,1,:),3).*xyz.sampleRate/10;
fet.data(:,2) = dot(uvec,nvec(:,2,:),3).*xyz.sampleRate/10;

% DIAGNOSTIC figure
% $$$ figure,
% $$$ ind = Trial.stc{'n'};
% $$$ ind = Trial.stc{'w'};
% $$$ hist2(fet(ind,:),...
% $$$       linspace(-100,100,50),...
% $$$       linspace(-100,100,50));

featureTitles = {};
featureDesc = {};
if nargout>1,
end


%---------------------------------------------------------------------------------------------------





