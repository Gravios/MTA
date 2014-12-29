classdef MTAModel
% MTAModel(model_base,flag) - Data structure holding marker names their relationships
%
%  CAUTION - All markers found in a *.vsk file must match those
%  found in  ~/data/config/MTA/MTAMarkers.mat. Add the marker name
%  manually or to MTAConfiguration.m (then run it) if it is not recognized.
%
%  model_base - string: name of the model (dependent upon flag)
%
%  flag - string: denotes whether the model should be constructed
%                 from the vicon *.vsk (vicon skelleton) file or created
%                 internally based on pre-configured model names  
%
%    accepted values: 
%
%      '-vsk' - the *.vsk file should be located in the xyz
%               directory associated with its vicon subsession
%
%                                |  session   |subsession|   model_base    |       
%               (e.g. ~/data/xyz/jg05-20120310/cof       /jg05-20120310-cof.vsk
%
%      '-mar' - select model based on the number of markers in the model
%               note: this flag creates a less informative model,
%                     always try to use the '-vsk' flag
%
%  Examples:
%    
%    create vsk model
%      Model = MTAModel(jg05-20120310-cof,'-vsk');
%
%    create mar model
%      Model = MTAModel(8,'-mar')
%

    properties (SetAccess = public)
        
        %name - string: Name of the model
        name

        %N - int: Number of markers in the model
        N

        %Connections - cellArray<-MTAStick: Array of MTAStick's containing connectivity information
        Connections

        %Markers - cellArray<-MTAMarker: Array of MTAMarker's containing marker information        
        Markers

        %imdMean - numericArray<-double: ignore this
        imdMean

        %imdStd - numericArray<-double: ignore this
        imdStd

        %index - int: xyz index where subject is in an ideal pose
        index
    end

    methods
        function Model = MTAModel(varargin)
            if isempty (varargin),
                Model.name = '';
                Model.Markers = {};
                Model.Connections = {};
                Model.imdMean = [];
                Model.imdStd  = [];
                Model.index = 0;
                return
            else
                [model_base,mflag] = DefaultArgs(varargin,{[],[]},true);
            end

            switch mflag
              case '-vsk'
                Markers = {};
                model_xml = parseXML(model_base);
                markers = xmlget(model_xml,'KinematicModel','MarkerSet','Markers',{'Marker','Attributes','NAME','RGB'});
                for i = 1:length(markers{1}),
                    Model.Markers{end+1} = MTAMarker(markers{1}{i}.Value,markers{2}{i}.Value);
                end

                Model.N = length(Model.Markers);

                numheadmar = 0;
                for i=1:Model.N,
                    if ~isempty(regexpi('head_front','^head'))
                        numheadmar = numheadmar+1;
                    end
                end
                
                Model.name = ['H' num2str(numheadmar) 'B' num2str(Model.N-numheadmar)];

                Model.Connections = {};
                sticks = xmlget(model_xml,'KinematicModel','MarkerSet','Sticks',{'Stick','Attributes','MARKER1','MARKER2','RGB'});
                for i = 1:length(sticks{1})
                    Model.Connections{end+1} = MTAStick(sticks{1}{i}.Value,sticks{2}{i}.Value,str2double(sticks{3}{i}.Value));
                end

                Model.imdMean = [];
                Model.imdStd  = [];
                Model.index = 0;
                return

              case '-mar'
                Model.N = length(model_base);
                Model.Markers = {};
                for i = 1:Model.N,
                    Model.Markers{end+1} = MTAMarker(model_base{i},'[0,0,0]');
                end

                load('MTAMarkerConnections.mat');

                switch Model.N
                  case 4
                    Model.name = 'H4';
                  case 5
                    Model.name = 'H5';
                  case 6
                    Model.name = 'H4B2';
                  case 7
                    Model.name = 'H5B2';
                  case 8
                    Model.name = 'H4B4';
                  case 9
                    Model.name = 'H5B4';
                  case 10
                    Model.name = 'H4B6';
                  case 11
                    Model.name = 'H5B6';
                  otherwise
                    error('model not registered: check MTAModel.m')
                end
                Model.Connections = {};
                Model.imdMean = [];
                Model.imdStd  = [];
                Model.index = 0;

                for j = 1:length(MTAMarkerConnections),
                    if sum(sign(Model.gmi(MTAMarkerConnections{j})))==2,
                        Model.Connections{end+1} = MTAStick(MTAMarkerConnections{j}{1},MTAMarkerConnections{j}{2},[0,0,0]);
                    end
                end
              
            end
        end 

        function markerIndex = gmi(Model,targetMarkers)
        %find marker index corresponding to marker(s) name
        %
        %Example:
        %  single marker index
        %    markerIndex = Session.Model.gmi('head_front')
        %
        %  multiple markers indecies 
        %    markerIndex = Session.Model.gmi('head_front')
        %
        %Output:
        %  markerIndex - numericArray<-int/int: index of markers
        %
            if ischar(targetMarkers)
                targetMarkers = {targetMarkers};
            end
            markerIndex = zeros(1,length(targetMarkers));
            for i = 1:length(targetMarkers)
                for j = 1:length(Model.Markers)
                    if strcmp(Model.Markers{j}.name,targetMarkers{i}),
                        markerIndex(i) = j;
                    end
                end
            end
        end

        function markerList = ml(Model)
        %List model marker names
        %
        %Example:
        %  markerList = Session.Model.ml();
        %
        %Output:
        %  markerList - cellArray<-string: list of marker names
        %
            markerList = {};
            for i = 1:Model.N,
                markerList{end+1} = Model.Markers{i}.name;
            end
        end

        function rigidBody = rb(Model,markerSet)
        %rigidBody = rb(Model,markerSet)
        %make a new model based on a subset of markers from an
        %existing model
            rigidBody = Model;
            rigidBody.N = length(markerSet);
            rigidBody.Markers = {};
            for i = 1:rigidBody.N,
                rigidBody.Markers{end+1} = Model.Markers{Model.gmi(markerSet{i})};
            end
            rigidBody.Connections = {};
            for i = 1:length(Model.Connections),
                if rigidBody.gmi(Model.Connections{i}.marker1)&&rigidBody.gmi(Model.Connections{i}.marker2),
                    rigidBody.Connections{end+1} = Model.Connections{i};
                end
            end
            
            try
                rigidBody.imdMean = [];
                rigidBody.imdStd = [];
                rigidBody.imdMean = Model.imdMean(Model.gmi(markerSet),Model.gmi(markerSet),:);
                rigidBody.imdStd = Model.imdStd(Model.gmi(markerSet),Model.gmi(markerSet),:);
            end
        end

        function DataCopy = copy(Data)
        % Make a copy of a handle object.
        % Instantiate new object of the same class.
            DataCopy = feval(class(Data));
            % Copy all non-hidden properties.
            p = properties(Data);
            for i = 1:length(p)
                if isa(Data.(p{i}),'MTAData'),
                    DataCopy.(p{i}) = Data.(p{i}).copy;
                else
                    DataCopy.(p{i}) = Data.(p{i});
                end
            end
        end

    

    end


    

    
    methods(Static)
        
        function validateMarkers(Model)
        %validateMarkers(Model)
        %check if all markers correspond to those saved in MTAMarkers.mat
        %
            valid = 0;
            load('MTAMarkers.mat');
            if ~exist('MTAMarkers','var'),
                error('configuration error: missing global variable valid_markers')
            end
            for i = 1:Model.N,
                for j = 1: length(MTAMarkers),
                    if strcmp(Model.Markers{i}.name,MTAMarkers{j}),
                        valid = valid + 1;
                        break
                    end
                end
            end
            if valid~=length(Model.markers)
                error('InvalidMarkerError: One or more markers do not exist in MTAMarkers.mat')
            end
        end
    
    end


end

