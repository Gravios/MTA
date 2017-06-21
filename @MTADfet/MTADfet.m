classdef MTADfet < MTAData
%MTADfet(path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key)
%
% <capspam>  THIS IS ALL WRONG - DO NOT READ FURTHER THAN THIS LINE
%
%  YOU ARE REALLY BAD AT FOLLOWING DIRECTONS </capspam> 
%
%  MTADxyz is a subclass of MTAData, which stores position data from a vicon
%  recording system.
%
%  NOTE: Future versions will allow other systems
%
%  Current Data Type: TimeSeries
%
%  Indexing (TimeSeries):
%    first dimension:    time, ':', numeric array of indicies or
%                              start and stop periods in an nx2 matrix 
%
%    second dimension:   marker, ':', numeric array of indicies or
%                              string corresponding to one of the model labels
%    
%    third dimension:    cartesian coordinates in R3 (x,y,z)
%
%    Indexing Example:
%       xy coordinates of 2 markers for all time
%       xy_head = xyz(:,{'head_back','head_front'},[1,2]);
%
%       Selected periods for z coordinates of a marker
%       z_head = xyz([1,300;400,1000],'head_front',3);
%
%
%  See also MTAData
    
    properties
        model
        data
    end
    
    methods
        function Data = MTADfet(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','fet',[],[],[]});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);   
            Data.model = model;
        end
        function Data = create(Data,varargin)
        %Data = create(Data,varargin)
        %not implemented in this version
        end
            

        function Data = embed(Data,win,overlap)
        %Data = embed(Data,win,overlap)
        %not implemented in this version
        end
        function center_of_mass = com(Data,Model,varargin)
        %center_of_mass = com(Session,Model)
        %
        %Model - MTAModel: MTA object holding marker information
        %
        %Examples:
        %  Find the center of mass of the session model
        %    center_of_mass = Session.com(Session.xyz.model);
        %
        %  Select a model based on a subset of markers from a larger model
        %    Model = Session.model.rb({'head_back','head_left','head_front','head_right'});
        %    center_of_mass = Session.com(Model);
        %   
        %
            [dim] = DefaultArgs(varargin,{[1:3]});
            center_of_mass = mean(Data.subsref(substruct('()',{':',Model.ml(),dim})),2);
        end

        
        function v = vel(Data,varargin)
        %v = vel(Data,varargin)
        %calculate the speed of marker(s)
        %[markers,dims] = DefaultArgs(varargin,{[1:Data.model.N],[1:size(Data.xyz,3)]});
        %zeros are added to the the beginning to for padding; reinterpolate
        %later to match xyz vector.
            [markers,dims] = DefaultArgs(varargin,{1:Data.model.N, 1:Data.size(3)});

            if iscell(markers)||ischar(markers),
                markers = Data.model.gmi(markers);
            elseif islogical(markers)
                markers = find(markers);
            end

            if ischar(markers),markers ={markers};end
            v = MTADxyz('data',cat(1,zeros([1,numel(markers)]),sqrt(sum(diff(Data.subsref(substruct('()',{':',markers,dims})),1,1).^2,3)).*Data.sampleRate./10),...
                        'sampleRate',Data.sampleRate);
            v.model = Data.model.copy;
            v.model.N = numel(markers);
            %v.model.Connections = 
            v.model.Markers = v.model.Markers(markers);
        end            

        function a = acc(Data,varargin)
        %a = acc(Data,varargin)
        %calculate the acceleration of marker(s)
        %[marker,dim] = DefaultArgs(varargin,{[1:Data.model.N],[1:size(Data.xyz,3)]});
            [markers,dims,padded] = DefaultArgs(varargin,{1:Data.model.N, 1:Data.size(3),1});
            vel = Data.vel(markers,dims);
            a = diff(vel.data,1);
            if padded==1
                a = MTADxyz('data',cat(1,a,a(end,:)),'sampleRate',Data.sampleRate);
            end
            a.model = vel.model.copy;
        end

         function Data = addMarker(Data,name,color,sticks,data)
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
             
             Marker = MTAMarker(name,color);
             Data.model.Markers{end+1} = Marker;
             Data.model.N = Data.model.N + 1;
             for i = 1:length(sticks),
                 Data.model.Connections{end+1} = MTAStick(sticks{i}{1},sticks{i}{2},sticks{i}{3});
             end
             Data.data = cat(2,Data.data,data);
         end


    end
    
    methods (Static)
        function Data = encapsulate(Trial,data,sampleRate,name,label,key)
        %function Data = encapsulate(Trial,data,sampleRate,name,label,key)
        %
        %  Argin:
        %    
        %    Trial: MTASession/MTATrial, motion tracking session or trial object
        %
        %    data:  Numeric,             feature data NxD 
        %                                N := number of samples
        %                                D := number of features
        %
        %    sampleRate: Numeric,        final sampleRate of output 
        %
        %    name:  String,              Name of the feature set
        %
        %    label: String,              Label/Tag for reference. (short and descriptive)
        %
        %    key:   Char,                Single character for keyboard shortcuts
        %
        %
            Data = MTADfet(Trial.spath,...
                           [],...
                           data,...
                           sampleRate,...
                           Trial.sync.copy,...
                           Trial.sync.data(1),...
                           [],'TimeSeries',[],...
                           name,label,key);                  
        end
        

    end
    
    
end