function Stc = label_bhv_homebase(Stc,Trial,varargin);


% DEFARGS -----------------------------------------------------------------------------------------
defargs = struct('stcMode',                     'msnn_ppsvd'                                   ...
);
[stcMode] = DefaultArgs(varargin,defargs,'--struct');
% -------------------------------------------------------------------------------------------------


if isempty(Stc),
    Stc = Trial.load('stc',stcMode);
end

state = Stc{'loc'};

[location,occupancy] = identify_homebase_locations(Trial);

xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQI');

xyz.addMarker('homebase',...     Name
              [1,0,0],...  Color
              {{'homebase','hcom',[0,0,255]}},... 
              permute(repmat([location(1,:),20],[size(xyz,1),1]),[1,3,2])...
);

fxyz = xyz.copy();
fxyz.filter('ButFilter',3,0.8,'low');


% COMPUTE inter marker angles and relative to homebase locations
ang = create(MTADang,Trial,fxyz);
angularOffsetHeadxHombase = MTADxyz('data',      circ_dist(ang(:,'acom','bcom',1),ang(:,'acom','homebase',1)),...
                                    'sampleRate',xyz.sampleRate);



nind = nniz(ang(:,'acom','homebase',3));
fdist = MTADxyz('data',      nan([size(ang,1),1]),...
                'sampleRate',xyz.sampleRate);
fdist.data(nind) = [diff(ButFilter(ang(nind,'acom','homebase',3),3,1/ang.sampleRate.*0.5,'low'));0];


homebase = state.copy;
homebase.clear;
homebase.label = 'locToHome';
homebase.key   = 'i';
homebase.updateFilename(Trial);

explore = state.copy;
explore.clear;
explore.label = 'locFromHome';
explore.key   = 'o';
explore.updateFilename(Trial);

for period = state.data',
    if sum(fdist(period'))>0,
        explore.data(end+1,:) = period';
    else
        homebase.data(end+1,:) = period';
    end
end

Stc.states{end+1} = homebase;
Stc.states{end+1} = explore;





%Stc.save(1);

