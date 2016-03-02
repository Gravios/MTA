function obj = MTAB_create_3d_object(Data,parent,handles,idx)

    stick_options = struct('width',      1,...
                           'erase',      'normal',... %{none|xor|background}
                           'line_style', '-');
    
    marker_options = struct('style','o',...
                            'size', 6 ,...
                            'erase', 'normal'); %{none|xor|background}
    
    % create button in MLobjectPane                   
    %guidata(hObject, handles);
    
    
    nSticks  = length(Data.model.Connections);
    nMarkers = Data.model.N;
    markerConnections=[];
    for i = 1:nSticks,
        markerConnections= cat(1,markerConnections,....
                               [Data.model.gmi(Data.model.Connections{i}.marker1),...
                                Data.model.gmi(Data.model.Connections{i}.marker2)]...
                              );
    end
    
    obj = struct('xyzpos',     Data.data,...
                 'nSticks',    nSticks,...                 
                 'sticks',     {cellfun(@feval,repmat({'animatedline'},[1,nSticks ]), 'UniformOutput',false)},...
                 'nMarkers',   nMarkers,...
                 'markers',    {cellfun(@feval,repmat({'animatedline'},[1,nMarkers]), 'UniformOutput',false)},...
                 'markerConnections', markerConnections);
    %% Translate connection map to pairs
    
    % SET sticks' properties
    for l=1:length(Data.model.Connections),
        obj.sticks{l}.Parent = parent;
        obj.sticks{l}.Tag = [Data.model.Connections{l}.marker1 ':' Data.model.Connections{l}.marker2];
        obj.sticks{l}.addpoints([obj.xyzpos(idx,markerConnections(l,1),1),obj.xyzpos(idx,markerConnections(l,2),1)],...
                                [obj.xyzpos(idx,markerConnections(l,1),2),obj.xyzpos(idx,markerConnections(l,2),2)],...
                                [obj.xyzpos(idx,markerConnections(l,1),3),obj.xyzpos(idx,markerConnections(l,2),3)]);
        obj.sticks{l}.LineWidth =  stick_options.width;
        %obj.sticks{l}.Color    , Session.model.Connections{l}.color./255, ...
        obj.sticks{l}.Visible = 'on';
        guidata(obj.sticks{l},handles);
    end
    
    % SET markers' properties
    for l=1:nMarkers,
        obj.markers{l}.Parent = parent;
        obj.markers{l}.Tag = Data.model.Markers{l}.name;
        obj.markers{l}.addpoints([obj.xyzpos(idx,l,1),obj.xyzpos(idx,l,1)],...
                                 [obj.xyzpos(idx,l,2),obj.xyzpos(idx,l,2)],...
                                 [obj.xyzpos(idx,l,3),obj.xyzpos(idx,l,3)]);
        obj.markers{l}.Marker           = marker_options.style;
        obj.markers{l}.MarkerEdgeColor  = Data.model.Markers{l}.color./255;
        obj.markers{l}.MarkerSize       = marker_options.size;
        obj.markers{l}.MarkerFaceColor  = Data.model.Markers{l}.color./255;
        obj.markers{l}.Visible          = 'on';
        guidata(obj.markers{l},handles);
    end
    