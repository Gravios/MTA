% req20170127.m 
% segment locomotive trajectories into in/out bound from homebase locations


Trial = MTATrial('jg05-20120309');

state = Trial.stc{'loc'};

[location,occpancy] = identify_homebase_locations(Trial);

xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD');

xyz.addMarker('homebase',...     Name
              [1,0,0],...  Color
              {{'homebase','hcom',[0,0,255]}},... 
              permute(repmat([location,0],[size(xyz,1),1]),[1,3,2])...
);



ang = create(MTADang,Trial,xyz);

angularOffsetHeadxHombase = circ_dist(ang(:,'acom','bcom',1),ang(:,'acom','homebase',1));


% $$$ figure,plot(angularOffsetHeadxHombase)
% $$$ figure,plot(ang(:,'acom','homebase',3))

nind = nniz(ang(:,'acom','homebase',3));
fdist = MTADxyz('data',nan([size(ang,1),1]),'sampleRate',xyz.sampleRate);
fdist.data(nind) = [diff(ButFilter(ang(nind,'acom','homebase',3),3,1/ang.sampleRate.*0.5,'low'));0];

% $$$ hold on,
% $$$ plot(fdist.data);
% $$$ 
% $$$ figure,
% $$$ plot(diff(fdist());

homebase = state.copy;
homebase.clear;
homebase.label = 'locToHome';
homebase.key   = 'i';
homebase.updateFilename(Trial);

explore = state.copy;
explore.clear;
explore.label = 'locFromHome';
explore.key   = 'o'
explore.updateFilename(Trial);

for period = state.data',
    if sum(fdist(period'))>0,
        explore.data(end+1,:) = period';
    else
        homebase.data(end+1,:) = period';
    end
end



figure,
eds = linspace(-pi,pi,100);
hm = bar(eds,histc(ang(homebase,'head_back','head_front',2),eds),'histc');
hm.FaceColor = 'c';
hm.FaceAlpha = 0.4;
hm.EdgeColor = 'c';
hm.EdgeAlpha = 0.4;
hold on,
hb = bar(eds,histc(ang(explore,'head_back','head_front',2),eds),'histc');
hb.FaceColor = 'r';
hb.FaceAlpha = 0.4;
hb.EdgeColor = 'r';
hb.EdgeAlpha = 0.4;


pfs_2d_states(Trial,[],{'loc','locToHome','locFromHome','rear','pause'}