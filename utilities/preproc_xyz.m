function [xyz,ss] = preproc_xyz(Trial,varargin)
%function [xyz,ss] = preproc_xyz(Trial,varargin)
%[procOpt] = DefaultArgs(varargin,{{}},true);
%
% An attempt to normalize marker positions along the spine
% between subjects using a 

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
      
      case 'spline_spine'
        ss = fet_spline_spine(Trial,'3dss',xyz);
        
        xyz.data(:,1:4,:) = ss(:,[5,35,65,95],:);
      
      case 'spline_spine_eqd'
        ss = fet_spline_spine(Trial,xyz);
        
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

      case 'SPLINE_SPINE_HEAD_EQD'

        if ~isempty(listFiles(Trial.name,'\.seh\.h')) && ~overwrite
            xyz = Trial.load('xyz','seh');
        else
            
            %Trial = MTASession.validate(Trial.filebase);
            xyz = Trial.load('xyz','trb');
            ss = fet_spline_spine(Trial,'3dssh',xyz,'overwrite',true);                

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
        ss = fet_spline_spine(Trial,'3dssh',xyz);                
      
      case 'SPLINE_SPINE_HEAD_EQD_NO_TRB'

        if ~isempty(listFiles(Trial.name,'\.seh\.h')) && ~overwrite
            xyz = Trial.load('xyz','seh');
        else
            
            %Trial = MTASession.validate(Trial.filebase);
            xyz = Trial.load('xyz');
            ss = fet_spline_spine(Trial,'3dssh',xyz,'overwrite',true);                

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
        ss = fet_spline_spine(Trial,'3dssh',xyz);                
       
      case 'load_trb_xyz'
        xyz = Trial.load('xyz','trb');
    
    end
    procOpts(1) = [];
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

              
              