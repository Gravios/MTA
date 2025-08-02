% FIGURES --------------------------------------------------------------------------------------------------------------
%
% What you will learn how to ... 
%   open a figure window.
%   populate a figure with axes
%   choose the units for the figure layout
%   some useful tricks to make the axes layout nice ... even with colorbars
%
% 1. Single Axes object
% 2. Many Axes objects
% 3. Mov'n Axes objects
% 4. Figure format 
% 5. Pesky Colorbar
%
% Protip: graphical handles inherit the hgsetget superclass, which allows you to use the get and set functions.
%         All objects are pass by reference.
%

%% 1. Single Axes object ---------------------------------------------------------------------------------------------------
% introducing the functions:
%     figure();
%     get();
%     axes();
%     gca();
%

hfig = figure(); % creates window for your plots

get(hfig)        % returns a struct with all of the properties and values of the figure

hax = axes();    % creates a new axes object within your figure

gca()            % get_current_axes -> returns the current Axes object

hfig.CurrentAxes % is the more direct version of gca();



%% 2. Many Axes objects ----------------------------------------------------------------------------------------------------
% introducing the functions:
%     clf();
%     subplot();
%     subplot2();
%     text();
%     num2str();

% UNCIVILIZED --------------------------------------------------------
clf();           % clears the current figure window

subplot(121);    % Create first of two axes in a horizontal arangement
subplot(122);    % Create second of two axes in a horizontal arangement


clf();           % clears the current figure window

subplot(211);    % Create first of two axes in a vertical arangement
subplot(212);    % Create second of two axes in a vertical arangement

clf();           % clears the current figure window

opts.text =                               ... keep blocky code stylish ... even in matlab
    {                                     ... 
                   'FontSize', 24,        ...
                 'FontWeight', 'bold',    ...
        'HorizontalAlignment', 'center'   ... 
    };

%                  xpos ypos text ...extra stuff
subplot(331); text( 0.5, 0.5, '1', opts.text{:});   % Create first of grid
subplot(332); text( 0.5, 0.5, '2', opts.text{:});   % Create second of grid
subplot(333); text( 0.5, 0.5, '3', opts.text{:});   % .
subplot(334); text( 0.5, 0.5, '4', opts.text{:});   % ... ah you get the idea
for s = 5:9,
    subplot(3,3,s);
    text( 0.5, 0.5, num2str(s), opts.text{:});
end
% NOTE the order ... do it!



% CIVILIZED ----------------------------------------------------------
clf();           % clears the current figure window
%        ny,nx, y, x
subplot2( 2, 1, 1, 1);    % Create first of two axes in a vertical arangement
subplot2( 2, 1, 2, 1);    % Create second of two axes in a vertical arangement



% ULTRA CIVILIZED ----------------------------------------------------
clf();
% FIRST row
%        ny,nx, y, x
subplot2( 3, 3, 1, 1);
subplot2( 3, 3, 1, 2);
subplot2( 3, 3, 1, 3);
% SECOND row
subplot2( 3, 3, 2, [1,2]);
subplot2( 3, 3, 2, 3);    
% THIRD row
subplot2( 3, 3, 3, [1,2,3]);


%% 3. Mov'n Axes objects ----------------------------------------------------------------------------------------------------
%
clf(); % clears the current figure window


opts.axes =                                   ...
    {                                         ...
            'Units', 'normalized',            ...
        'LineWidth', 1                        ...
    };

opts.axesPositions =                          ...
    { ...  x    y    w    h  
        [0.2, 0.6, 0.2, 0.2],                 ... top left
        [0.6, 0.6, 0.2, 0.2],                 ... top right
        [0.2, 0.3, 0.6, 0.1]                  ... bottom
    };

        
hax = gobjects([3,1]);
%        ny,nx, y, x
for a = 1:numel(opts.axesPositions)
    hax(a) = axes('Position', opts.axesPositions{a}, ...
                  opts.axes{:});
end






%% 4. Figure format (my version) ---------------------------------------------------------------------------------------
% 
hfig = figure();
[hfig, opts, fax, sax] = set_figure_layout(hfig,'A4','portrait');

%opts.subplot.width





%% 5. Pesky Colorbar ---------------------------------------------------------------------------------------------------

imdata = imread('ngc6543a.jpg'); % read in some image data

nColors = size(imdata,3);
nx = nColors;
ny = 2;

hfig = figure();

for chan = 1:nColors,
    subplot2( 2, 3, 1, chan);
    imagesc(imdata(:,:,chan));
    colorbar();

    subplot2( 2, 3, 2, chan);    % Create second of two axes in a vertical arangement
    histogram(                   ... of a single color channel from a slice of imdata which was reshapde into a column vector
        reshape(imdata(:,:,chan),...
                [],              ...
                1),              ...
        'BinEdges',1:255);
end

% The color bar breaks the visual flow sometimes, so you have to put it in its place!

clf();           % clears the current figure window

for chan = 1:nColors,
    hax = subplot2( 2, 3, 1, chan); % store Axes object handle
    imagesc(imdata(:,:,chan));
    cax = colorbar(hax);
    cax.Position(1) = sum(hax.Position([1,3])); % 

    subplot2( 2, 3, 2, chan);    % Create second of two axes in a vertical arangement
    histogram(                   ... of a single color channel from a slice of imdata which was reshapde into a column vector
        reshape(imdata(:,:,chan),...
                [],              ...
                1),              ...
        'BinEdges',1:255);
end
% That looks better



