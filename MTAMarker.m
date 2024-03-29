classdef MTAMarker
% MTAMarker
% Class which holds general information and graphical properties of
% the markers.
%
% The naming convention of markers is restricted by the C3D format
% to a set length of ##. The name of a marker should be a two word
% or a double acronym separated by an underscore. Choose words
% whos initials are unique (necessary for clarity when the short
% function is called).
%
    properties (SetAccess = public)
        name
        color
    end
    methods 
        function Marker = MTAMarker(name,color)
            Marker.name = name;
            if ischar(color),
                Marker.color = str2num(color);
                if any(Marker.color>1)
                    Marker.color = Marker.color/255;
                end
            else
                Marker.color = color(:)';
            end
        end
        
        function strout = short(Marker)
            ts = strjoin(regexpi(Marker.name,'(^.)|(_.)','match'),'');
            strout = upper(ts([1:2:numel(ts)]));
        end
        
    end
end


