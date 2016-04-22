function [xyz,ss] = preproc_xyz(Trial,varargin)
%function [xyz,ss] = preproc_xyz(Trial,varargin)
%[procOpt] = DefaultArgs(varargin,{{}},true);

[procOpts] = DefaultArgs(varargin,{{}},true);


%Trial = MTATrial.validate(Trial);


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

if ~iscell(procOpts), procOpts = {procOpts}; end
    
while ~isempty(procOpts)
    switch procOpts{1}
      
      case 'spline_spine'
        try 
            ss = Trial.load('fet','3dss');
            ss.resample(xyz);
        catch
            txyz = xyz.data;
            pnts = zeros([xyz.size(1),105,3]);
            for ind = 1:xyz.size(1),
                try
                    pnts(ind,:,:) = fnplt(cscvn(sq(txyz(ind,:,:))'))';
                end
            end
            name = '3d spline interpolated spine'; label = '3dss'; key = 's';
            ss = MTADfet.encapsulate(Trial,...
                                     pnts,...
                                     xyz.sampleRate,...
                                     name,label,key);
            ss.updateFilename(Trial);
            ss.save;
        end
        xyz.data(:,1:4,:) = ss(:,[5,35,65,95],:);
      
      case 'spline_spine_eqd'
        try 
            ss = Trial.load('fet','3dss');
            ss.resample(xyz);
        catch
            txyz = xyz.data;
            pnts = zeros([xyz.size(1),105,3]);
            for ind = 1:xyz.size(1),
                try
                    pnts(ind,:,:) = fnplt(cscvn(sq(txyz(ind,1:4,:))'))';
                end
            end
            name = '3d spline interpolated spine'; label = '3dss'; key = 's';
            ss = MTADfet.encapsulate(Trial,...
                                     pnts,...
                                     xyz.sampleRate,...
                                     name,label,key);
            ss.updateFilename(Trial);
            ss.save;
        end
         
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
