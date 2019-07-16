function markerPoints = label_frames(varargin)
% function markerPoints = label_frames(varargin)
%    
% Zooming and Panning are functional
%
% INPUTS
%
%     path:        String,    location of images
%     markers:     CellArray, names of the markers
%     markerPoint: Matrix,    < numFrames, numMarkers, 2 > - in case you need to take a break
%       
% Instructions
%    
%    left click: set marker     
%    right click: remove marker 
%    middle click: next frame
%    keyboard:                  
%        n : next marker        
%        p : previous marker    
%        q : quit               
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('path',                pwd(),                                                   ...
                 'markers',             {{'tail_dist','tail_prox','spine_lower',                 ...
                                          'pelvis_root','spine_middle','spine_upper',            ...
                                          'head_back','head_left','head_right'}},                ...
                 'markerPoints',        [],                                                      ...
                 'autofocusFlag',       true                                                     ...
);
[path,markers,markerPoints,autofocusFlag] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------    


% $$$ path = pwd();
% $$$ path = '/storage2/nogay/labeled-data/Camera13029';

files = dir(path);
files(1:2) = [];

hfig = figure(2039320);
clf(hfig);
hfig.Units = 'centimeters';
hfig.Position = [5,5,20,12];

% SET image axes
hax = axes('Position',[0.1,0.1,0.6,0.8]);
hax.XTick = [];
hax.YTick = [];
hold(hax,'on');
pax = pan();
zax = zoom();
axis(hax,'ij');
axis(hax,'tight');

% SET background axes
fax = axes('Units','Normalized','Position',[0,0,1,1],'Visible','off');
xlim(fax,[0,1]);
ylim(fax,[0,1]);

% PRINT help section
text(0.75,0.35,{'left click: set marker     '});
text(0.75,0.3, {'right click: remove marker '});
text(0.75,0.25,{'middle click: next frame   '});
text(0.75,0.2 ,{'keyboard:                  '});
text(0.75,0.15,{'    n : next marker        '});
text(0.75,0.1, {'    p : previous marker    '});
text(0.75,0.05,{'    q : quit               '});

% CREATE marker table
tab = uicontrol('style','listbox',...
                'String',markers,...
                'Max',1,...
                'units','normalized',...
                'position',[0.7,0.4,0.3,0.5],...
                'background','w',...
                'Value',1);


numMarkers = numel(markers);

cax = gobjects([1,numel(markers)]);

currentMarkerInd = 1;
if isempty(markerPoints)
    markerPoints = zeros([numel(files),numMarkers,2]);
end

img = [];
for f = 1:numel(files),
    imgtmp = imread(fullfile(path,files(f).name),'png');
    img(:,:,f) = imgtmp(:,:,1);
end
bkg = median(img,3);


quitFlag = false;
for f = find(any(any(markerPoints,3),2),1,'last'):numel(files),
    cla(hax);
    axes(hax);
    imagesc(hax,img(:,:,f));
    if autofocusFlag
        rat = imgaussfilt(img(:,:,f)-bkg,18);
        [~,ratPos] = max(rat(:));        
        [y,x] = ind2sub(size(rat),ratPos);
        xlim([x-200,x+200]);
        ylim([y-200,y+200]);        
    end
    exitFlag = false;    
    hfig.CurrentCharacter = 'x';    
    currentMarkerInd = 1;    
    while ~exitFlag,    
        waitforbuttonpress();
        if ~strcmp(hfig.CurrentCharacter,'x'),
            switch hfig.CurrentCharacter
              case 'n', % next marker
                currentMarkerInd = currentMarkerInd+1;                
                if numMarkers < currentMarkerInd,
                    currentMarkerInd = 1;
                end
                tab.Value = currentMarkerInd;
                
              case 'p', % previous marker
                currentMarkerInd = currentMarkerInd-1;                
                if 1 > currentMarkerInd,
                    currentMarkerInd = numMarkers;
                end
                tab.Value = currentMarkerInd;
                
              case 'f', % quit
                exitFlag = true;
                
              case 'q'
                exitFlag = true;
                quitFlag = true;
                
            end
            continue
        end        
            
        if strcmp(zax.Enable,'on') | strcmp(pax.Enable,'on'),
            continue
        end
        
        selectionType = get(hfig, 'selectiontype');
        switch selectionType
          case 'normal', % set new marker
            markerPoint = hax.CurrentPoint;
            markerPoints(f,currentMarkerInd,:) = round(markerPoint(1,1:2));
            currentMarkerInd = currentMarkerInd+1;                
            if numMarkers < currentMarkerInd,
                currentMarkerInd = 1;
            end
            tab.Value = currentMarkerInd;
            
          case 'alt', % delete current marker
            currentMarkerInd = currentMarkerInd-1;                
            if 1 > currentMarkerInd,
                currentMarkerInd = numMarkers;
            end
            tab.Value = currentMarkerInd;
            markerPoints(f,currentMarkerInd,:) = zeros([1,1,2]);
          case 'extend',
            exitFlag = true;            
        end
        
        % RESET current character be
        hfig.CurrentCharacter = 'x';
        
        % UPDATE circles
        for c = 1:numMarkers
            delete(cax(c));
            cax(c) = circle(markerPoints(f,c,1),markerPoints(f,c,2),10,'m');
        end
    end
    if quitFlag,
        break
    end
end


