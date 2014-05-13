function out = draw_lines(h,mode)
% Click the left mouse button to define a point
% Drag the mouse to draw a line to the next point and
% left click again
% Right click the mouse to stop drawing

% 
set(h,'WindowButtonDownFcn',@wbdcb)
ah = gca;
setappdata(h,'pntCnt',0);
setappdata(h,'lnCnt',0);
setappdata(h,'mode',mode);
setappdata(h,'state',true);
exlines = findobj('-regexp','tag','dll');  if ~isempty(exlines),delete(exlines),end
expoints = findobj('-regexp','tag','dlp'); if ~isempty(expoints),delete(expoints),end

   function wbdcb(src,evnt)
      if strcmp(get(src,'SelectionType'),'normal')       
         [x,y,str] = disp_point(ah);

         lnCnt = getappdata(src,'lnCnt');
         hl = line('XData',x,'YData',y,'Marker','.');
         set(hl,'tag',['dll',num2str(lnCnt+1)])
         setappdata(src,'lntCnt',lnCnt+1);

         pntCnt = getappdata(src,'pntCnt');
         text(x,y,str,'VerticalAlignment','bottom','tag',['dlp',num2str(pntCnt+1)]);drawnow
         setappdata(src,'pntCnt',pntCnt+1);

         set(src,'WindowButtonMotionFcn',@wbmcb)
      elseif strcmp(get(src,'SelectionType'),'alt')
         set(src,'WindowButtonMotionFcn','')
         [x,y,str] = disp_point(ah);
         pntCnt = getappdata(src,'pntCnt');
         text(x,y,str,'VerticalAlignment','bottom','tag',['dlp',num2str(pntCnt+1)]);drawnow
         setappdata(src,'pntCnt',pntCnt+1);
         
         switch getappdata(src,'mode')
           case 'line_fit'
             hl = findobj('-regexp','tag','dll');
             for i = 1:numel(hl)
                 prams(i,:) = polyfit(get(hl(i),'XData'),get(hl(i),'YData'),1);
             end
             setappdata(src,'out',prams);
             setappdata(h,'state',false);
         end

      end
      function wbmcb(src,evnt)
         [xn,yn,str] = disp_point(ah);
         xdat = [x,xn];
         ydat = [y,yn];
         set(hl,'XData',xdat,'YData',ydat);
      end  
   end
   function [x,y,str] = disp_point(ah)
      cp = get(ah,'CurrentPoint');  
      x = cp(1,1);y = cp(1,2);
      str = ['(',num2str(x,'%0.3g'),', ',num2str(y,'%0.3g'),')'];    
   end

while getappdata(h,'state'),pause(.1);end

out = getappdata(h,'out');

end