function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_href_HXYZ_nose(Trial,varargin)
% function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_href_HXY(Trial,varargin)
% 
% varargin:
%     newSampleRate: numeric,  (Trial.xyz.sampleRate) - sample rate of xyz data
%     normalize:     logical,  (false)                - covert each feature to z-score
%     procOpts:      CellARY,  ({'SPLINE_SPINE_HEAD_EQD'}), - preprocessing options
%     filterCutoff:  numeric,  (4) - lowpass filter in Hz
%     theta:         numeric,  (0) - angle to rotate head coordinates
%     markers:       CellARY,  {'head_back','head_left','head_front','head_right'} 

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
fxyz = preproc_xyz(Trial,procOpts);
fxyz.resample(newSampleRate);
fxyz.data(~nniz(fxyz),:,:) = 0;
fxyz.filter('ButFilter',4,filterCutoff,'low');

xyz = preproc_xyz(Trial,procOpts);
xyz.resample(newSampleRate);
xyz.data(~nniz(xyz),:,:) = 0;
xyz.filter('ButFilter',4,20,'low');


hcom = fxyz(:,'hcom',[1,2,3]);

% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(fxyz(:,'head_back',:)-hcom,fxyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 
nm = nz+hcom;
fxyz.addMarker('htx',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
ny = cross(fxyz(:,'htx',:)-hcom,fxyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny+hcom;
fxyz.addMarker('hrx',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
nx = cross(fxyz(:,'hrx',:)-hcom,fxyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx+hcom;
fxyz.addMarker('hbx',  [0.5,1,0.5],[],nm);


nb = cat(2,nx,ny,nz);
% $$$ nb = cat(2,nx,ny,nz,hcom);
% $$$ nb = cat(3,nb,zeros([size(nb,1),size(nb,2),1]));
% $$$ nb(:,4,4) = 1;
% $$$ nb = permute(nb,[1,3,2]);
%sq(nb(1010,:,:))*inv(sq(nb(1000,:,:)))*sq(nb(1010,:,:))

%hcom = xyz(:,'hcom',[1,2,3]);

% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 
nm = nz+hcom;
xyz.addMarker('htx',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny+hcom;
xyz.addMarker('hrx',  [0.5,1,0.5],[],nm);

% GENERATE orthogonal basis, origin: head's center of mass
nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx+hcom;
xyz.addMarker('hbx',  [0.5,1,0.5],[],nm);

newx = multiprod(nb,sq(nx),[2,3],2);
newy = multiprod(nb,sq(ny),[2,3],2);
nbh = cat(2,nx,ny,nz);
% $$$ nbh = cat(2,nx,ny,nz,hcom);
% $$$ nbh = cat(3,nbh,zeros([size(nb,1),size(nb,2),1]));
% $$$ nb(:,4,4) = 1;
% $$$ nb = permute(nb,[1,3,2]);


figure,
subplot(211);
hold('on');
plot(xts,newx(:,2))
plot(xts,newy(:,1))
% $$$ plot(xts,newx(:,2))
% $$$ plot(xts,newx(:,3))
colormap('jet');
subplot(212); 
    hold('on');
    plotSTC(stc,1,'text',{'theta','sit','groom','pause','turn','walk','rear'},'kymcgbr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');




% GENERATE orthogonal basis, origin: head's center of mass
nvec = xyz(:,'nose',[1,2,3])-hcom;
nvec = bsxfun(@rdivide,nvec,sqrt(sum((nvec).^2,3))); 

nvec = cat(2,nvec,cross(nvec(:,1,:),nvec(:,2,:),3));


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





