function plot_distrbs(Sessions,var1,var2,varargin)
[name,states,v1nbins,v2nbins,stc_mode] = DefaultArgs(varargin,{'jpdf_headHeight_bSpeed','rw',30,30,'auto_wbhr')

% $$$ %% Test Vars
% $$$ %                  SessionName   TrialName DataServer
% $$$ Sessions = {{'er01-20110719',     'all', 'bach'},... CA3
% $$$                {'er01-20110721',     'all', 'bach'},... CA3
% $$$                {'er06-20130612',     'all',  'cin'},... CA1
% $$$                {'er06-20130613', 'all-cof',  'cin'},... CA1
% $$$                {'er06-20130614', 'all-cof',  'cin'},... CA1
% $$$                {'jg04-20120129',     'all', 'bach'},... CA1
% $$$                {'jg04-20120130',     'all', 'bach'},... CA1
% $$$                {'jg05-20120309',     'all', 'bach'},... CA1
% $$$                {'jg05-20120310',     'all', 'bach'},... CA1
% $$$                {'jg05-20120311',     'all', 'bach'},... CA1
% $$$                {'jg05-20120315',     'all', 'bach'},... CA1
% $$$                {'jg05-20120317',     'all', 'bach'},... CA2???
% $$$                {'jg05-20120324',     'all', 'bach'}}; % CA3
% $$$ stc_mode = 'auto_wbhr';
% $$$ %s = 12;
% $$$ %% end - Test vars


nses = numel(Sessions);

for s = 1:nses

    MTAstartup('cin',Sessions{s}{3})
    Trial = MTATrial(Sessions{s}{1},Sessions{s}{2});
    Trial.stc.updateMode(stc_mode);
    try 
        Trial.stc.load;
    catch
        if Trial.stc.isempty,
            try
                ds = load(Trial.stc.fpath);
                Trial.stc.states = ds.states;
                Trial.stc.save(1);
            catch
                Trial = labelBhv(Trial,Trial.stc);
            end
        end
    end
    %[state,type] = DefaultArgs(varargin,{[],'linear'})

% $$$ %% Test Vars
% $$$ states = 'rw';
% $$$ v1nbins = 30;
% $$$ v2nbins = 30;
% $$$ xyz = Trial.xyz.copy;
% $$$ xyz.load(Trial);
% $$$ xyz.filter(gtwin(.5,xyz.sampleRate));
% $$$ vel = MTADxyz('data',log10([0;Trial.vel('spine_lower',[1,2])]),'sampleRate',xyz.sampleRate);
% $$$ %% end - Test vars
% $$$ var1 = vel;
% $$$ var2 = MTADxyz('data',log10(xyz(:,7,3)),'sampleRate',xyz.sampleRate);

no_zeros = var1~=0&var2~=0&var1~=nan&var2~=nan&var1~=inf&var2~=inf&var1~=-inf&var2~=-inf;


%% Plot JPDF(v1,v2), distribution 
figh = figure(2323),
subplot(121),hold on
ind = no_zeros;

v1lim =prctile(var1(ind),[1,99]);
v2lim =prctile(var2(ind),[1,99]);
v1edgs = linspace(v1lim(1),v1lim(2),v1nbins);
v2edgs = linspace(v2lim(1),v2lim(2),v2nbins);

N = hist2([var1(ind),var2(ind)],v1edgs,v2edgs);
grd = cell(1,2);
[grd{:}] = ndgrid(mean([v1edgs(1:end-1)',v1edgs(2:end)'],2),mean([v2edgs(1:end-1)',v2edgs(2:end)'],2));
grd{3}=N/sum(N(find(N)));
contour(grd{:},20,'k');

%% Test Vars
ylabel('head height log(mm)')
xlabel('body speed log(mm/s)')
%% end - Test vars

subplot(122),hold on
c = 'rgbcmy';
cd = {'red','green','blue','cyan','magenta','yellow'};
slab = {};
for i = 1:numel(states)
    ind = Trial.stc{states(i)};
    slab{end+1} = ind.label; 
    N = hist2([var1(ind),var2(ind)],v1edgs,v2edgs);
    grd{3}=N/sum(N(find(N)));
    contour(grd{:},20,c(i));
end

%% Test Vars
ylabel('head height log(mm)')
xlabel('body speed log(mm/s)')
%% end - Test vars

legend(slab{:},'Location','southwest')
suptitle(['Contour JPDF ' Trial.filebase ])

reportfig('/gpfs01/sirota/bach/homes/gravio/figures',figh,'jpdfc',[],Trial.filebase,100)
end