function [cpnts]=ClusterPP(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[figureId] = DefaultArgs(varargin,{gcf});

setappdata(figureId,'oldFigName',get(figureId,'Name'));
set(figureId,'Name',mfilename)


tax = gca;
setappdata(figureId,'target_axes', tax);

h = findobj(tax,'Type','line');
x = get(h,'Xdata');
y = get(h,'Ydata');
if isempty(x) || isempty(y)
   warning(['The current figure is empty or does not exist'])
   cpnts=[];
   close(figureId)
   return
else
    setappdata(figureId,'x',x); clear('x');
    setappdata(figureId,'y',y); clear('y');
    setappdata(figureId,'cluster_points',zeros(size(getappdata(figureId,'x'))));
end

setappdata(figureId,'axis_hold_status',ishold(tax));
oldAxisTitle = copyobj(get(tax,'Title'),tax);
set(oldAxisTitle,'Visible','off');
setappdata(figureId,'axis_title',oldAxisTitle);
hold(tax,'on');

setappdata(figureId,'local_x',[]);
setappdata(figureId,'local_y',[]);
setappdata(figureId,'temp_x', []);
setappdata(figureId,'temp_y', []);
setappdata(figureId,'sh',     []);
setappdata(figureId,'clu_sh', []);
setappdata(figureId,'area_sh',[]);


setappdata(figureId,'fcluster',true);
setappdata(figureId,'cborders',[]  );
setappdata(figureId,'clustercoordinates',[]);
setappdata(figureId,'n_cluster'         , 1);
setappdata(figureId,'non_empty_cluster' , 1);
setappdata(figureId,'exit'              , false);

set(figureId, 'KeyPressFcn',          @KeyPressFun);
set(figureId, 'WindowButtonDownFcn'  ,@clicker1);
set(figureId, 'WindowButtonMotionFcn',[]);



uiwait(figureId);

exit = false;


while ~exit,
% $$$     if non_empty_cluster==0
% $$$         display(['Last selection was empty!'])
% $$$         display(['Select the group of points for the cluster ' num2str(n_cluster)])
% $$$     else
% $$$         display(['Select the group of points for the cluster ' num2str(n_cluster)])
% $$$     end
% $$$     if n_cluster==1
% $$$         display(['Please select the limits of the cluster. Left mouse click confirms a vertex. Rigth mouse click confirms the area cover by the vertices.'])
% $$$         display(['The title of the figure indicates the current coordenates of the pointer'])
% $$$     end
% $$$ 
% $$$     display(['Left mouse click to set the first vertex'])


    data = getappdata(figureId);
    clustercoordinates=[data.local_x,data.local_y];
    set(getappdata(figureId,'sh'),'XData',[data.local_x,data.temp_x],'YData',[data.local_y,data.temp_y]);
    drawnow;
    set(figureId, 'WindowButtonDownFcn', @clicker1);
    set(figureId, 'WindowButtonMotionFcn',@mouseMove1);
 
    uiwait(figureId);
 
    data = getappdata(figureId);
    exit = data.exit;
    

end


data = getappdata(figureId);



figure;
line(data.x(data.cluster_points==0),data.y(data.cluster_points==0),'LineStyle','.','Color',[.8 .8 .8])
hold on
color=colormap(hsv(max(data.cluster_points)));
for cl=1:max(data.cluster_points),
    line(data.x(data.cluster_points==cl),data.y(data.cluster_points==cl),...
         'LineStyle',  '.',                                    ...
         'Color',      color(cl,:))

end

cpnts = data.cluster_points;
set(figureId, 'KeyPressFcn',[]);
end



function mouseMove1 (figureId,eventdata)

C = get(getappdata(figureId,'target_axes'), 'CurrentPoint');
title(getappdata(figureId,'target_axes'), ['(X,Y) = (', num2str(C(1,1)), ', ',num2str(C(1,2)), ')']);
setappdata(figureId,'temp_x',C(1,1));
setappdata(figureId,'temp_y',C(1,2));

if ~isempty(getappdata(figureId,'sh')),
    set(getappdata(figureId,'sh'),...
        'XData',[getappdata(figureId,'local_x'), getappdata(figureId,'temp_x')],...
        'YData',[getappdata(figureId,'local_y'), getappdata(figureId,'temp_y')]);
    drawnow;
else
    setappdata(figureId,'temp_x', []);
    setappdata(figureId,'temp_y', []);
% $$$     setappdata(figureId,'sh',...
% $$$         'XData',[data.local_x, data.temp_x],...
% $$$         'YData',[data.local_y, data.temp_y]);

end

end


function clicker1 (figureId,eventdata)
data = getappdata(figureId);

S=get(figureId, 'selectiontype');

switch(S)
    case 'normal' % 'normal' for left moue button   
        C = get (data.target_axes, 'CurrentPoint');            
        setappdata(figureId,'temp_x',C(1,1));
        setappdata(figureId,'temp_y',C(1,2));
        setappdata(figureId,'clustercoordinates',...
                            [getappdata(figureId,'clustercoordinates');...
                            getappdata(figureId,'temp_x'),getappdata(figureId,'temp_y')]);

        data = getappdata(figureId);
        

        if isempty(data.sh),

            setappdata(figureId,'sh',line([data.temp_x, data.temp_x],...
                                          [data.temp_y, data.temp_y],...
                                          'LineStyle'  ,'-',    ...
                                          'Color'      ,[0,0,1],...
                                          'linewidth'  ,2));
        else
            
            set(data.sh,'XData',[data.local_x, data.temp_x],...
                        'YData',[data.local_y, data.temp_y]);
        end

        drawnow;

        setappdata(figureId,'area_sh',[data.area_sh copyobj(getappdata(figureId,'sh'),...
                                                            get(getappdata(figureId,'sh'),'parent'))]);
        setappdata(figureId,'local_x',data.temp_x);
        setappdata(figureId,'local_y',data.temp_y);

        
    case 'alt'  % 'alt' for right mouse button

        C = get (gca, 'CurrentPoint');
        setappdata(figureId,'temp_x',C(1,1));
        setappdata(figureId,'temp_y',C(1,2));
        setappdata(figureId,'clustercoordinates',...
                            [getappdata(figureId,'clustercoordinates');...
                            getappdata(figureId,'temp_x'),getappdata(figureId,'temp_y')]);

        data = getappdata(figureId);


        set(data.sh,'XData',[data.local_x, data.temp_x],...
                    'YData',[data.local_y, data.temp_y]);


        setappdata(figureId,'area_sh',[data.area_sh copyobj(getappdata(figureId,'sh'),...
                                                            get(getappdata(figureId,'sh'),'parent'))]);
        
        setappdata(figureId,'local_x',data.temp_x);
        setappdata(figureId,'local_y',data.temp_y);


        cl=inpolygon(data.x,data.y,data.clustercoordinates(:,1),data.clustercoordinates(:,2));
        ch=plot(data.x(cl),data.y(cl),'xk');
        setappdata(figureId,'clu_sh',[data.clu_sh ch]);

        set(figureId, 'WindowButtonMotionFcn', []);


        ash = getappdata(figureId,'area_sh');
        for i = 1:numel(ash),delete(ash(i));end
        setappdata(figureId,'area_sh',[]);

        setappdata(figureId,'clustercoordinates',[]);
        if nnz(cl)
            data.cluster_points(cl)=data.n_cluster;
            setappdata(figureId,'cluster_points',data.cluster_points);
            setappdata(figureId,'non_empty_cluster',1);
        else
            warning(['ATTENTION: empty cluster'])
            setappdata(figureId,'non_empty_cluster',0);
        end
                
        delete(getappdata(figureId,'sh'));
        setappdata(figureId,'sh',[]);

        drawnow;

        setappdata(figureId,'n_cluster',data.n_cluster+data.non_empty_cluster);

%     case 'extended' % 'extend' for middle mouse button
%         
%     case 'open' % 'open' on double click

        
  otherwise 
    %eventdata.Key
        set (figureId, 'WindowButtonMotionFcn', []);
        set(figureId, 'WindowButtonDownFcn', []);
        
        delete(data.clu_sh)
        return


        
end
        
uiresume(figureId);        


end

function KeyPressFun (figureId,eventdata)
    key = get(figureId,'CurrentCharacter');
    if ~isempty(key),
        switch double(key)
          case 27 % '^['
            setappdata(figureId,'exit',true);
            delete(getappdata(figureId,'sh'));
            setappdata(figureId,'sh',[]);
            set(figureId, 'KeyPressFcn',          []);
            set(figureId, 'WindowButtonDownFcn'  ,[]);
            set(figureId, 'WindowButtonMotionFcn',[]);
            delete(getappdata(figureId,'clu_sh'));
            
            tax = getappdata(figureId,'target_axes');
            ahold = getappdata(figureId,'axis_hold_status');
            if ~ahold, hold(tax,'off'); end                
            delete(get(tax,'Title'));
            oldAxisTitle = getappdata(figureId,'axis_title');
            set(tax,'Title',oldAxisTitle);
            set(oldAxisTitle,'Visible','on');
            set(figureId,'Name',getappdata(figureId,'oldFigName'));
        end
    end
    
    uiresume(figureId);        
end
