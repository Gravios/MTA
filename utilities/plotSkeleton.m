function [markers,sticks,markerConnections,traj] = plotSkeleton(Trial,xyz,ind,varargin)
%function [markers,sticks,markerConnections] = plotSkeleton(xyz,ind,varargin)
%[plot_type,ang,t_ante,t_post] = DefaultArgs(varargin,{'line',[],0,0},true);

default_args = ...
{... plot_type
    'line',...
 ... ang
     [],...
 ... t_per
    [0,0],...
 ... t_makrers
    {'spine_lower','pelvis_root','spine_middle','spine_upper','head_left','head_right'},...
 ... hax 
    gca ...
};


[plot_type,ang,t_per,t_markers,hax] = DefaultArgs(varargin,default_args);


if strcmp(class(hax),'matlab.ui.Figure'),hax = gca; end

Rx = @(theta) [1,0,0;0,cos(theta),-sin(theta);0,sin(theta),cos(theta)];
Ry = @(theta) [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
Rz = @(theta) [cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];

traj = [];
markerConnections=[];
for i = 1:length(xyz.model.Connections)
    markerConnections= cat(1,markerConnections,....
                           [xyz.model.gmi(xyz.model.Connections{i}.marker1),xyz.model.gmi(xyz.model.Connections{i}.marker2)]);
end

sticks = repmat({0},1,length(xyz.model.Connections));

for l=1:length(xyz.model.Connections), 
    switch plot_type
      case 'line'
        sticks{l} =animatedline([xyz.data(ind(end),markerConnections(l,1),1),...
                                    xyz.data(ind(end),markerConnections(l,2),1)],...
                                [xyz.data(ind(end),markerConnections(l,1),2),...
                                    xyz.data(ind(end),markerConnections(l,2),2)],...
                                [xyz.data(ind(end),markerConnections(l,1),3),...
                                    xyz.data(ind(end),markerConnections(l,2),3)]);
       set(sticks{l},'Parent',hax,...
                     'Tag',[xyz.model.Connections{l}.marker1 ':' xyz.model.Connections{l}.marker2],...
                     'LineWidth', 2, ... stick_options.width,
                  ...%'Color'    , xyz.model.Connections{l}.color./255, ...
                     'Visible','on',...
                     'Parent',hax);


         % dotted trajectory before and/or after the ind
         if any(t_per~=0), 
             for i= 1:numel(t_markers)
                 traj(end+1)=animatedline(xyz(ind+(t_per(1):1:t_per(2)),t_markers{i},1),...
                                          xyz(ind+(t_per(1):1:t_per(2)),t_markers{i},2),...
                                          xyz(ind+(t_per(1):1:t_per(2)),t_markers{i},3));
                 set(traj(end),'LineStyle',':')
                 set(traj(end),'Marker',   '.',...
                               'color',[.5,.5,.5],...
                               'tag',t_markers{i});
                 set(traj(end),'MarkerSize',4)
                 set(traj(end),'Parent',hax);
             end
         end

         
      case 'surface'
        %compute rotatios
        mang = sq(ang(ind(end),markerConnections(l,1),markerConnections(l,2),:));
        mxyz = permute(mean(xyz(ind(end),markerConnections(l,:),:),2),[1,3,2]);
        cxyz = [];
        [cxyz(:,:,1),cxyz(:,:,2),cxyz(:,:,3)] = cylinder(2,50);
        cxyz(:,:,3) = (cxyz(:,:,3)-.5).*mang(3);
        cxyz = reshape(cxyz,[],3);
        cxyz = cxyz*(Ry(pi/2-mang(2))*Rz(pi-mang(1))*Rx(0));
        cxyz = bsxfun(@plus,cxyz,mxyz);
        cxyz = reshape(cxyz,2,[],3);

        sticks{l} = surf(cxyz(:,:,1),cxyz(:,:,2),cxyz(:,:,3),...
                         'facecolor',[0,0,1],... %xyz.model.Connections{l}.color./255,
                         'edgecolor','none','Parent',hax);        

    end

end


%% Markers
marker_options = {};
marker_options.style ='o' ;
marker_options.size = 8 ;
marker_options.erase = 'none'; %{none|xor|background}
markers = repmat({0},1,xyz.model.N);
for l=1:xyz.model.N,
    switch plot_type
      case 'line'
        markers{l} =animatedline([xyz.data(ind(end),l,1),xyz.data(ind(end),l,1)],...
                         [xyz.data(ind(end),l,2),xyz.data(ind(end),l,2)],...
                         [xyz.data(ind(end),l,3),xyz.data(ind(end),l,3)]);
        set(markers{l},'Parent',hax, ...
                       'Tag',xyz.model.Markers{l}.name,...
                       'Marker'         , marker_options.style, ...
                       'MarkerEdgeColor', xyz.model.Markers{l}.color./255, ...
                       'MarkerSize'     , marker_options.size, ...
                       'MarkerFaceColor', xyz.model.Markers{l}.color./255, ...
                       'Visible','on');

      case 'surface'

        [xx,yy,zz] = meshgrid(-15:15,-15:15,-15:15);
        rr = sqrt(xx.^2 + yy.^2 + zz.^2);
        xx = xx+xyz.data(ind(end),l,1);
        yy = yy+xyz.data(ind(end),l,2);
        zz = zz+xyz.data(ind(end),l,3);
        %# calculate distance from center of the cube

        %# create the isosurface by thresholding at a iso-value of 10
        markers{l} = patch(isosurface(xx,yy,zz,rr,10));
        isonormals(xx,yy,zz,rr,markers{l});
        %set(markers{l}, 'VertexNormals', isonormals(xx,yy,zz,rr,markers{l}));
        set(markers{l}, 'FaceColor', xyz.model.Markers{l}.color./255,...
                        'EdgeColor', 'none','Parent',hax);
% $$$                         'XData',get(markers{l},'XData')+xyz.data(ind(end),l,1),...
% $$$                         'YData',get(markers{l},'YData')+xyz.data(ind(end),l,2),...
% $$$                         'ZData',get(markers{l},'ZData')+xyz.data(ind(end),l,3));

    end
end

% $$$ bom = xyz.com(xyz.model);
% $$$ [xf,yf,zf] = meshgrid([-250:250]+bom(ind,1,1),[-250:250]+bom(ind,1,2),0);
% $$$ hfloor = surface(xf,yf,zf,repmat(permute([0,0,0],[1,3,2]),[size(xf),1]),'edgecolor','none');

daspect(hax,[1 1 1]);
c = camlight;
c.Parent = hax;
lighting(hax,'phong')
