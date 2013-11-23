classdef MTAMarker
    properties (SetAccess = public)
        name
        color
    end
    methods 
        function Marker = MTAMarker(name,color)
            Marker.name = name;
            Marker.color = str2num(color);
        end
    end
end

