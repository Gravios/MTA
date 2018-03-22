function [fet,featureTitles,featureDesc] = fet_template(Trial,varargin)
%function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis(Trial,varargin)
% varargin:
%    newSampleRate    numeric    (12)    -    sample rate of output:fet 
%    normalize        Logical    (false) -    apply unity transform
%    procOpts         CellStr    ({'SPLINE_SPINE_HEAD_EQD'})
%                                        -    preprocessing options for xyz
% varargout:
%    fet              MTADfet
%    featureTitles    CellStr
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', 12,                       ...
                 'normalize'    , false,                    ...
                 'procOpts'     , {{'SPLINE_SPINE_HEAD_EQD'}});

[newSampleRate,normalize,procOpts] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'req20160310_selected_features','fet_mis','m');                  

% XYZ preprocessed 
xyz = preproc_xyz(Trial,procOpts);
resample(xyz,newSampleRate);


% CAT feature
fet.data = [ ];

% SET zero where xyz is zero
fet.data(~nniz(xyz),:)=0;

% SET labels and descriptions 
featureTitles = {};
featureDesc = {};
if nargout>1,
    % 1.
    featureTitles(end+1) = {'feature_name'};    
    featureDesc(end+1)   = {['Description of feature']};
end


%---------------------------------------------------------------------------------------------------