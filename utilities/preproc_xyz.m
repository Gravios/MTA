function [xyz,ss] = preproc_xyz(Trial,varargin)
%function [xyz,ss] = preproc_xyz(Trial,varargin)
%[procOpts,overwrite] = DefaultArgs(varargin,{{},false},true);
%
% An attempt to normalize marker positions along the spine
% between subjects using a 
%
% Current Opts:
%    NONE
%    SPLINE_SPINE
%    SPLINE_SPINE_EQD
%    SPLINE_SPINE_HEAD_EQI
%    SPLINE_SPINE_HEAD_EQD
%    SPLINE_SPINE_HEAD_EQD_NO_TRB
%    LOAD_TRB_XYZ
% 

[procOpts,overwrite] = DefaultArgs(varargin,{{},false},true);
Trial = MTATrial.validate(Trial);

% DELBLOCK if version>=3,
if isempty(Trial.fet),
    Trial.fet = MTADfet(Trial.spath,...
                        [],...
                        [],...
                        [],...
                        Trial.sync.copy,...
                        Trial.sync.data(1),...
                        []);                  
end
% DELBLOCK end


% XYZ Positions of Markers
xyz = Trial.load('xyz');
ss = [];
%ss = fet_spline_spine(Trial,'3dss',xyz);

if ~iscell(procOpts), procOpts = {procOpts}; end
    
while ~isempty(procOpts)
    switch procOpts{1}
      
      case 'SPLINE_SPINE'
        ss = fet_spline_spine(Trial,'3dss','pos');
        
        xyz.data(:,1:4,:) = ss(:,[5,35,65,95],:);
      
      case 'SPLINE_SPINE_EQD'
        ss = fet_spline_spine(Trial,'pos');
        
        spineLength = MTADxyz('data',sqrt(sum(diff(ss.data,1,2).^2,3)),'sampleRate',xyz.sampleRate);
        totalSpineLength = sum( spineLength.data ,2);
        cumSpineLength = cumsum( spineLength.data ,2);
        meanSpineLength = nanmean(totalSpineLength);        
        nNewMarkers = 4;
        xs = zeros([spineLength.size(1),4]);
        for t = 1:spineLength.size(1)-1,
            [~,xi,~] = NearestNeighbour(cumSpineLength(t,:),...
                                          [0:nNewMarkers-1].*totalSpineLength(t)/(nNewMarkers-1),'both');
            xyz.data(t,1:4,:) = ss(t,xi,:);
            xs(t,:) = xi;
        end
        xyz.label = 'sed';
        xyz.key  = 'e';
        xyz.name = 'spline_spine_eqd';
        xyz.updateFilename(Trial);      
        xyz.save

      case 'SPLINE_SPINE_HEAD_EQI'
        
        targetMarkers = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
        numMarkers = numel(targetMarkers);

        if ~isempty(list_files(Trial.name,'.sehs.h')) && ~overwrite
            xyz = Trial.load('xyz','sehs');
        else
            
            
            xyz = Trial.load('xyz','trb');
            %xyz.filter('ButFilter',3,50,'low');
            ss = fet_spline_spine(Trial,'3dssh','trb');                

            % COM head Center of Mass
            xyz.addMarker('hcom',...     Name
                          [.7,0,.7],...  Color
                          {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
                           {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
                           {'head_front','hcom',[0,0,255]},...
                           {'head_right','hcom',[0,0,255]}},... 
                             xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
            
            spineSegmentLength = MTADxyz('data',sqrt(sum(diff(ss.data,1,2).^2,3)),'sampleRate',xyz.sampleRate);
            nind = find(nniz(xyz));
            ssn = ss.copy();
            ssn.data = zeros([size(ss,1),size(ss,2)-6,size(ss,3)]);
            markerInd = zeros([size(ss,1),numMarkers-2]);
            for t = nind'
                nzind = [true,spineSegmentLength(t,:)]>1e-10;
                zind = find(~nzind)-1;
                nzind([zind]) = ~nzind([zind]);    
                ssn.data(t,:,:) = ss(t,nzind,:);
                markerInd(t,:) = zind-[1:numel(zind)].*2+spineSegmentLength(t,[zind-1])./sum(reshape(spineSegmentLength(t,[zind-1,zind+1]),3,2),2)';
            end
            %spineSegmentLengthNew = MTADxyz('data',sqrt(sum(diff(ssn.data,1,2).^2,3)),'sampleRate',xyz.sampleRate);
            baseInd = 100/(numMarkers-1);

            samplePoints = 1:size(ssn,2);
            for m = 1:numMarkers-2,
                medianMarkerIndOffset = baseInd*m-median(markerInd(nind,m));
                xMarInd = xyz.model.gmi(targetMarkers{m+1});
                for d = 1:size(ssn,3),
                    for t = nind',
                        xyz.data(t,xMarInd,d) = interp1(samplePoints,ssn(t,:,d),markerInd(t,m)...
                                                        +medianMarkerIndOffset,'spline');
                    end
                end
            end

            if isa(Trial,'MTATrial')
                xyz.sync = Trial.sync.copy;
                xyz.sync.sync = Trial.sync.copy;
            end
            xyz.label = 'sehs';
            xyz.key  = 'h';
            xyz.name = 'spline_spine_head_eqi_smooth';
            xyz.updateFilename(Trial);      
            xyz.save;
            Trial = MTATrial.validate(Trial.filebase);
            xyz = Trial.load('xyz','sehs');
        end
        
        if nargout>1,            
            ss = fet_spline_spine(Trial,'3dssh','trb'); 
% $$$             nind = find(nniz(xyz));
% $$$             ssn = ss.copy();
% $$$             ssn.data = zeros([size(ss,1),size(ss,2)-6,size(ss,3)]);
% $$$             for t = nind'
% $$$                 nzind = [true,spineSegmentLength(t,:)]>1e-6;
% $$$                 zind = find(~nzind)-1;
% $$$                 nzind([zind]) = ~nzind([zind]);    
% $$$                 ssn.data(t,:,:) = ss(t,nzind,:);
% $$$             end
% $$$             ss = ssn;
        end

        
        
      case 'SPLINE_SPINE_HEAD_EQD'

        if ~isempty(list_files(Trial.name,'.seh.h')) && ~overwrite
            xyz = Trial.load('xyz','seh');
        else
            
            %Trial = MTASession.validate(Trial.filebase);
            xyz = Trial.load('xyz','trb');
            ss = fet_spline_spine(Trial,'3dssh','trb','overwrite',true);                

            % COM head Center of Mass
            xyz.addMarker('hcom',...     Name
                          [.7,0,.7],...  Color
                          {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
                           {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
                           {'head_front','hcom',[0,0,255]},...
                           {'head_right','hcom',[0,0,255]}},... 
                             xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
            
            mid = xyz.model.gmi({'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'});
            spineLength = MTADxyz('data',sqrt(sum(diff(ss.data,1,2).^2,3)),'sampleRate',xyz.sampleRate);
            totalSpineLength = sum( spineLength.data ,2);
            cumSpineLength = cumsum( spineLength.data ,2);
            meanSpineLength = mean(totalSpineLength(nniz(totalSpineLength)));
            nNewMarkers = 5;
            for t = 1:spineLength.size(1)-1,
                [~,xi,~] = NearestNeighbour(cumSpineLength(t,:),...
                                            [0:nNewMarkers-1].*totalSpineLength(t)/(nNewMarkers-1),'both');
                xyz.data(t,mid,:) = ss(t,xi,:);
            end
            if isa(Trial,'MTATrial')
                xyz.sync = Trial.sync.copy;
                xyz.sync.sync = Trial.sync.copy;
            end
            xyz.label = 'seh';
            xyz.key  = 'h';
            xyz.name = 'spline_spine_head_eqd';
            xyz.updateFilename(Trial);      
            xyz.save 
            Trial = MTATrial.validate(Trial.filebase);
            xyz = Trial.load('xyz','seh');
        end
        
      
      case 'SPLINE_SPINE_HEAD_EQD_NO_TRB'

        if ~isempty(list_files(Trial.name,'.seh.h')) && ~overwrite
            xyz = Trial.load('xyz','seh');
        else
            
            %Trial = MTASession.validate(Trial.filebase);
            xyz = Trial.load('xyz');
            ss = fet_spline_spine(Trial,'3dssh','trb','overwrite',true);                

            % COM head Center of Mass
            xyz.addMarker('hcom',...     Name
                          [.7,0,.7],...  Color
                          {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
                           {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
                           {'head_front','hcom',[0,0,255]},...
                           {'head_right','hcom',[0,0,255]}},... 
                             xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
            
            mid = xyz.model.gmi({'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'});
            spineLength = MTADxyz('data',sqrt(sum(diff(ss.data,1,2).^2,3)),'sampleRate',xyz.sampleRate);
            totalSpineLength = sum( spineLength.data ,2);
            cumSpineLength = cumsum( spineLength.data ,2);
            meanSpineLength = mean(totalSpineLength(nniz(totalSpineLength)));
            nNewMarkers = 5;
            for t = 1:spineLength.size(1)-1,
                [~,xi,~] = NearestNeighbour(cumSpineLength(t,:),...
                                            [0:nNewMarkers-1].*totalSpineLength(t)/(nNewMarkers-1),'both');
                xyz.data(t,mid,:) = ss(t,xi,:);
            end
            if isa(Trial,'MTATrial')
                xyz.sync = Trial.sync.copy;
                xyz.sync.sync = Trial.sync.copy;
            end
            xyz.label = 'seh';
            xyz.key  = 'h';
            xyz.name = 'spline_spine_head_eqd';
            xyz.updateFilename(Trial);      
            xyz.save 
            Trial = MTATrial.validate(Trial.filebase);
            xyz = Trial.load('xyz','seh');
        end
       
      case 'LOAD_TRB_XYZ'
        xyz = Trial.load('xyz','trb');
    
    end
    procOpts(1) = [];
end

if nargout>1,
    ss = fet_spline_spine(Trial,'3dssh','trb');                
end


% COM lower Body Center of Mass
xyz.addMarker('bcom',...     Name
              [.7,0,.7],...  Color
              {{'spine_lower', 'acom',[0,0,255]},... Sticks to visually connect
               {'pelvis_root', 'acom',[0,0,255]},... new marker to skeleton
               {'spine_middle','acom',[0,0,255]}},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));

% COM head Center of Mass
xyz.addMarker('hcom',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
               {'head_front','hcom',[0,0,255]},...
               {'head_right','hcom',[0,0,255]}},... 
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

% COM subject Center of Mass              
xyz.addMarker('acom',...    Name
              [.7,0,.7],... Color
              {{'spine_lower', 'acom',[0,0,255]},... Sticks to visually connect
               {'pelvis_root', 'acom',[0,0,255]},... new marker to skeleton
               {'spine_middle','acom',[0,0,255]},...
               {'spine_upper', 'acom',[0,0,255]},...
               {'head_back',   'acom',[0,0,255]},...
               {'head_front',  'acom',[0,0,255]}},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front'})));

              
              