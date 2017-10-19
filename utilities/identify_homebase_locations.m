function [location,occupancy,normalizeOccupancyMap,xyBins] = identify_homebase_locations(Trial,varargin)
%function homebaseLocation = identify_homebase_locations(Trial)
% Returns putative homebase locations of rodent subject
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('state',                   'pause+sit+groom',                                   ...
                 'stc',                     'msnn_ppsvd',                                        ...
                 'type',                    'xy',                                                ...
                 'binSize',                 20,                                                  ...
                 'smoothingWeights',        [5,5],                                               ...
                 'xyzProcOpts',             {{'SPLINE_SPINE_HEAD_EQI'}}                          ...
);
[state,stc,type,binSize,smoothingWeights,xyzProcOpts] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN --------------------------------------------------------------------------------------------- 

stc = Trial.load('stc');
xyz = preproc_xyz(Trial,xyzProcOpts);
state = stc{state};
state.data(diff(state.data,1,2)<state.sampleRate*10,:) = [];
state.data = [state.data(:,1),state.data(:,1)+round(state.sampleRate)];

[smoothedOccupancyMap,xyBins] = generate_occupancy_map(Trial,xyz,state,type,binSize,smoothingWeights,[],0.001);

 normalizeOccupancyMap = smoothedOccupancyMap{1}./nansum(smoothedOccupancyMap{1}(:));

[location,occupancy] = LocalMinimaN(-normalizeOccupancyMap,...
                                    -max(normalizeOccupancyMap(:))/2,...
                                    200/diff(xyBins{1}(1:2)));
location(:,1) = xyBins{1}(location(:,1));
location(:,2) = xyBins{2}(location(:,2));

% END MAIN -----------------------------------------------------------------------------------------