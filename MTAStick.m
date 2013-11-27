classdef MTAStick
% MTAStick(marker1,marker2,color) - Data structure to visualize connections/spatial organizations of markers
%
% example:   MTAStick('head_front','head_back',[0,0,1])
%

    properties
        %marker1 - string: marker descriptor (e.g.'bodypart_position')
        marker1

        %marker2 - string: marker descriptor (e.g.'bodypart_position')
        marker2

        %color - numericArray: rgb values range - [0,1] (e.g.[0.4,1,0])
        color
    end

    methods
        function Stick = MTAStick(marker1,marker2,color)
            Stick.marker1 = marker1;
            Stick.marker2 = marker2;
            Stick.color = color;
        end
    end
end