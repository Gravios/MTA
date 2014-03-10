classdef MTAMarker
    properties (SetAccess = public)
        name
        color
    end
    methods 
        function Marker = MTAMarker(name,color)
            Marker.name = name;
            if ischar(color),
                Marker.color = str2double(color);
            end
        end
    end
end

