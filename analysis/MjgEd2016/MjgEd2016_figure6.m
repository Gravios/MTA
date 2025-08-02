% MjgEd2016_figure6
%
% Behavioral composition of exploration


% $$$ OwnDir = '/storage/gravio/nextcloud/';
% $$$ FigDir = 'Shared/Behavior Paper/Figures/Figure_6/parts';
% $$$ mkdir(fullfile(OwnDir,FigDir));
% $$$ Trials   = af(@(t)  MTATrial.validate(t),  get_session_list('hand_labeled'));
% $$$ stc      = cf(@(t)  t.load('stc'),         Trials);
% $$$ xyz      = cf(@(t)  preproc_xyz(t,'SPLINE_SPINE_HEAD_EQI'),  Trials);
% $$$ stcm      = cf(@(t)  t.load('stc','msnn_ppsvd'),         Trials);
% $$$ %for s = 1:numel(stcm), stcm{s}.states(9:end) = [];end
% $$$ stcm = cf(@(s)  reduce_stc_to_loc(s),  stcm);
% $$$ stcm = cf(@(s,t) label_bhv_homebase(s,t), stcm, Trials);
% $$$ label_bhv_homebase(stcm{2},Trials{2})


Trial = MTATrial.validate('jg04-20120128.cof.all');
stc = Trial.load('stc','msnn_ppsvd');
%stc = reduce_stc_to_loc(stc);
%label_bhv_homebase(stc,Trial);
xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQI');



stopLocations = [];
sts = 'pause+groom+sit';
sper = stc{s};
for p = sper.data',
    if diff(p)>1240,
        stopLocations(end+1,:) = sq(mean(xyz(p','acom',[1,2])));
    end
end



% FIND homebase locations
[location,occupancy,map,bins] = identify_homebase_locations(Trial,'smoothingWeights',[5,5]);

% ADD homebase locations to xyz object
for l = 1:size(location,1),
    xyz.addMarker(['homebase',num2str(l)],...     Name
                  [1,0,0],...                     Color
                  {{['homebase',num2str(l)],'hcom',[0,0,255]}},... 
                  permute(repmat([location(l,:),100],[size(xyz,1),1]),[1,3,2])...
    );
end


% FILTER lowpass < 0.05 Hz xyz
fxyz = xyz.copy();
fxyz.filter('ButFilter',3,0.05,'low');
% COMPUTE speed of lowpass filtered markers
fvxy = fxyz.vel([],[1,2]);
fvxy.data(fvxy.data<1e-5) = 1e-6;
fvxy.data = log10(fvxy.data);
fvxy.filter('ButFilter',3,0.05,'low');

% COMPUTE distance between rat and hombases
fxyz = xyz.copy();
fxyz.filter('ButFilter',3,0.1,'low');
ang = create(MTADang,Trial,xyz);


figure,histogram2(stopLocations(:,1),stopLocations(:,2),bins{1},bins{2},'DisplayStyle','tile')
colormap jet

figure();
subplot2(5,5,1:3,1:5);
hold('on');
plot(fvxy(:,1));
Lines([],-1,'k');
plot(ang(:,'homebase1','acom',3)./100);
Lines([],3,'k');
subplot2(5,5,4:5,1:5);
plotSTC(stc);
linkaxes(findobj(gcf,'Type','Axes'),'x');

figure();
hold('on');
imagesc(bins{1},bins{2},map');
axis('xy');
plot(stopLocations(:,1),stopLocations(:,2),'m.'); xlim([-500,500]);ylim([-500,500]);
scatter(location(:,1),location(:,2),50,[1,0,0]);
scatter(location(:,1),location(:,2),200,[1,0,0]);
scatter(location(:,1),location(:,2),500,[1,0,0]);
scatter(location(:,1),location(:,2),1000,[1,0,0]);
scatter(location(:,1),location(:,2),1750,[1,0,0]);

figure();
plot(ang(:,'homebase1','acom',3))

plot(xyz(homebase,'hcom',1),xyz(homebase,'hcom',2),'.b')
plot(xyz(explore,'hcom',1),xyz(explore,'hcom',2),'.r')



% $$$ figure();hold('on');
% $$$ 
% $$$ figure,
% $$$ subplot(211); hist(ang(homebase,'head_back','head_front',2),linspace(-1.5,1.5,100));
% $$$ subplot(212); hist(ang(explore,'head_back','head_front',2),linspace(-1.5,1.5,100));
% $$$ 
% $$$ vxy = fxyz.vel([],[1,2]);
% $$$ vxy.data(vxy.data<1e-3) = 1e-3;
% $$$ vxy.data = log10(vxy.data);
% $$$ 
% $$$ figure,
% $$$ subplot(211); histogram(vxy(homebase,'hcom'),linspace(-0.5,1.75,100));
% $$$ subplot(212); histogram(vxy(explore,'hcom'), linspace(-0.5,1.75,100));
