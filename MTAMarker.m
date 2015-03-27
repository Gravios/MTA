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
            end
        end
        
        function strout = short(Marker)
            t = 'head_front';
            p = '(^.)|(_.)';
            ts = strjoin(regexpi(t,p,'match'),'');
            strout = upper(ts([1:2:numel(ts)]));
        end
        
    end
end


