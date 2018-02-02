function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_href_H(Trial,varargin)
% function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis(Trial,varargin)
% 
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate'  , Trial.xyz.sampleRate,                                        ...
                 'normalize'      , false,                                                       ...
                 'procOpts'       , {''},                                                        ...
                 'rotationAxisInd', 2,                                                           ...
                 'theta'          , 0                                                            ...
);
[newSampleRate,normalize,procOpts,rotationAxisInd,theta] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              Trial.xyz.sampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'head referenced position and motion','fet_href_H','h');                  

% PREPROC xyz
% FILTER xyz

xyz = preproc_xyz(Trial,procOpts);
xyz.data(~nniz(xyz),:,:) = 0;
xyz.filter('RectFilter',3,4);
rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
fhcom = ButFilter(hcom,3,[3]./(xyz.sampleRate/2),'low');

% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 
ny = cross(nz,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nx = cross(ny,nz);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));

if theta ~= 0,
    evec = cat(2,nx,ny,nz);
    j =1:3;
    headNorm = bsxfun(@rdivide,sq(evec(:,rotationAxisInd,:)),sqrt(sum(evec(:,rotationAxisInd,:).^2,3)));
    headKron = reshape(repmat(headNorm',3,1).*headNorm(:,j(ones(3,1),:)).',[3,3,size(headNorm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    headCPM = reshape(headNorm(:,k)',3,3,size(headNorm,1)).*repmat(j,[1,1,size(headNorm,1)]);

% CREATE rotation matrix
    j = 1:3;
    headRotMat = cos(theta)*repmat(eye(3),[1,1,size(headNorm,1)])...
        +sin(theta)*headCPM...
        +(1-cos(theta))*headKron;

% SET matrix
    ovec = evec(:,1,:);
% ROTATE Basis
    nx = permute(sum(headRotMat.*permute(reshape(nx(:,j(ones(3,1),:)),[size(headNorm,1),3,3]),[2,3,1]),2),[3,2,1]);
    ny = permute(sum(headRotMat.*permute(reshape(ny(:,j(ones(3,1),:)),[size(headNorm,1),3,3]),[2,3,1]),2),[3,2,1]);
    nz = permute(sum(headRotMat.*permute(reshape(nz(:,j(ones(3,1),:)),[size(headNorm,1),3,3]),[2,3,1]),2),[3,2,1]);
end


% COMPUTE hcom projection onto head reference
uvec = circshift(hcom,-3)-circshift(hcom,3);
fet.data(:,1) = dot(uvec,nx,3);
fet.data(:,2) = dot(uvec,ny,3);
fet.data(:,3) = dot(uvec,nz,3);

fet.resample(newSampleRate);
featureTitles = {};
featureDesc = {};
if nargout>1,
end


%---------------------------------------------------------------------------------------------------





% $$$ 
% $$$ hfig = figure(20392010);
% $$$ set(hfig,'units','centimeters')
% $$$ set(hfig,'Position',[11 0 32 25])
% $$$ set(hfig,'PaperPositionMode','auto');
% $$$ clf();
% $$$ hold('on');
% $$$ sind = 98000;
% $$$ hline = gco;
% $$$ mcolors = 'bgrm';
% $$$ mlineWidth = ones([1,numel(mcolors)]); mlineWidth(end)=3;
% $$$ mlineMarker = repmat({'none'},[1,numel(mcolors)]);
% $$$ mlineMarker{end} = 'o';
% $$$ for m = 1:size(nvec,2),
% $$$     hline(end+1) = line([0,nvec(sind,m,1)],[0,nvec(sind,m,2)],[0,nvec(sind,m,3)],...
% $$$                         'Color',mcolors(m),...
% $$$                         'LineWidth',mlineWidth(m),...
% $$$                         'Marker',mlineMarker{m});
% $$$ end
% $$$ daspect([1,1,1]);

