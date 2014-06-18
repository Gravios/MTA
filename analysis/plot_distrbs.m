function plot_distrbs(Sessions)


%% Test Vars
%                  SessionName   TrialName DataServer
SessionList = {{'er01-20110719',     'all', 'bach'},... CA3
               {'er01-20110721',     'all', 'bach'},... CA3
               {'er06-20130612',     'all',  'cin'},... CA1
               {'er06-20130613', 'all-cof',  'cin'},... CA1
               {'er06-20130614', 'all-cof',  'cin'},... CA1
               {'jg04-20120129',     'all', 'bach'},... CA1
               {'jg04-20120130',     'all', 'bach'},... CA1
               {'jg05-20120309',     'all', 'bach'},... CA1
               {'jg05-20120310',     'all', 'bach'},... CA1
               {'jg05-20120311',     'all', 'bach'},... CA1
               {'jg05-20120315',     'all', 'bach'},... CA1
               {'jg05-20120317',     'all', 'bach'},... CA2???
               {'jg05-20120324',     'all', 'bach'}}; % CA3
stc_mode = 'auto_wbhr';
s = 12;
%% end - Test vars


nses = numel(SessionList);

for s = 1:nses

MTAstartup('cin',SessionList{s}{3})
Trial = MTATrial(SessionList{s}{1},SessionList{s}{2});
Trial.stc.updateMode(stc_mode);
Trial.stc.load;

%[state,type] = DefaultArgs(varargin,{[],'linear'})

%% Test Vars
states = 'wgl';
v1nbins = 30;
v2nbins = 30;
xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.5,xyz.sampleRate));
vel = MTADxyz('data',log10([0;Trial.vel('spine_lower',[1,2])]),'sampleRate',xyz.sampleRate);
%% end - Test vars


var1 = vel;
var2 = MTADxyz('data',log10(xyz(:,7,3)),'sampleRate',xyz.sampleRate);

no_zeros = var1~=0&var2~=0&var1~=nan&var2~=nan&var1~=inf&var2~=inf&var1~=-inf&var2~=-inf;

figure,
%% Plot JPDF(v1,v2), distribution 
ind = no_zeros;
v1lim =prctile(var1(ind),[2,98]);
v2lim =prctile(var2(ind),[2,98]);
v1edgs = linspace(v1lim(1),v1lim(2),v1nbins);
v2edgs = linspace(v2lim(1),v2lim(2),v2nbins);
N = hist2([var1(ind),var2(ind)],v1edgs,v2edgs);
%caxis([0,2000]);
contour(v1edgs,v2edgs,N);


for i = 1:numel(states)
figure
ind = Trial.stc{states(i)};
hist2([var1(ind),var2(ind)],v1edgs,v2edgs);colorbar
end

end