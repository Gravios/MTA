function Data = addMarker(Data,varargin)
%Session = addMarker(Session,name,color,sticks,xyz)
%modify Model and xyz by adding an additional marker
%
%  name - string: marker descriptor (e.g.'bodypart_position')
%  color - numericArray: rgb values range - [0,1] (e.g.[0.4,1,0])
%  sticks - cellArray: information used to visualize connections/spatial
%                      organizations of markers
%                      {'marker_name1','marker_name2',[stick_color]}
%                      (e.g. {'head_front','head_back',[0,0,1]})
%
%             marker_name1 - string: marker descriptor (e.g.'bodypart_position')
%             marker_name2 - string: marker descriptor (e.g.'bodypart_position')
%             stick_color - numericArray: rgb values range - [0,1] (e.g.[0.4,1,0])
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('name',   'newMarker',                                                          ...
                 'color',  [.7,1,.7],                                                            ...
                 'sticks', {{}},                                                                 ...
                 'data',   []                                                                    ...
);        
[name,color,sticks,data] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------  


%check if marker exists
markerIndex = Data.model.gmi(name);
if markerIndex,
    %warning('MTA:MTADxyz:addMarker:MarkerExists, Overwriting Marker');
    Data.data(:,markerIndex,:) = data;
else
    Marker = MTAMarker(name,color);
    Data.model.Markers{end+1} = Marker;
    Data.model.N = Data.model.N + 1;
    if ~isempty(sticks),
        for i = 1:length(sticks),
            Data.model.Connections{end+1} = MTAStick(sticks{i}{1},sticks{i}{2},sticks{i}{3});
        end
    end
    Data.data = cat(2,Data.data,data);
    
end


