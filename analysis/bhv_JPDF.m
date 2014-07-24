% Test Vars
%                  SessionName   TrialName DataServer
Sessions = {{'er01-20110719',     'all', 'bach'},... CA3
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

%s = 12;
%% end - Test vars

%% Test Vars
%                  SessionName   TrialName DataServer
%Sessions = {{'jg05-20120310',     'all', 'bach'}};
% $$$ Sessions = {{'jg05-20120324',     'all', 'bach'}};
%Sessions = {{'er06-20130613', 'all-cof',  'cin'}};
%Sessions = {{'er01-20110719',     'all', 'bach'}};

nses = numel(Sessions);

for s = 1:nses

    MTAstartup('cin',Sessions{s}{3},false)
    try
        Trial = MTATrial(Sessions{s}{1},Sessions{s}{2});
    catch 
        continue
    end
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
                %Trial.stc.states{end+1} = theta(Trial,'sts2epoch');
                %Trial.stc.save(1);
            end
        end
    end

% $$$ Trial.load('xyz');
% $$$ Trial.xyz.filter(gtwin(1.2,Trial.xyz.sampleRate));
% $$$ vel = Trial.xyz.vel(1,[1,2]);
% $$$ plot(vel.data);Lines(Trial.stc{'w'}(:),[],'k');
% $$$ title([Trial.name,' ',num2str(s)])
% $$$ pause(15)
% $$$ end


    %% Test Vars
    states = 'rw';
    v1nbins = 30;
    v2nbins = 30;
    xyz = Trial.xyz.copy;
    xyz.load(Trial);
    xyz.filter(gtwin(.5,xyz.sampleRate));
    vel = xyz.vel('spine_lower',[1,2]);
    vel.data = log10(vel.data);
    ang = Trial.ang.copy;
    ang.create(Trial,xyz);
    %% end - Test vars
    var1 = vel; xlab = 'body speed log(mm/s)';
    %var1 = MTADxyz('data',ang(:,5,7,2),'sampleRate',xyz.sampleRate); xlab = 'head pitch (rad)';
    var2 = MTADxyz('data',log10(xyz(:,7,3)),'sampleRate',xyz.sampleRate); ylab = 'head height log(mm)';         


    no_zeros = nniz(var1)&nniz(var2);


    %% Plot JPDF(v1,v2), distribution 
    figh = figure(2323),
    subplot(131),hold on
    ind = no_zeros;

    v1lim =prctile(var1(ind),[1,99]);
    v2lim =prctile(var2(ind),[1,99]);
    v1edgs = linspace(v1lim(1),v1lim(2),v1nbins);
    v2edgs = linspace(v2lim(1),v2lim(2),v2nbins);


    N = hist2([var1(ind),var2(ind)],v1edgs,v2edgs);
    imagesc(v1edgs,v2edgs,N',[0,mean(N(nniz(N(:))))*3])
    xlim(v1edgs([1,end]))
    ylim(v2edgs([1,end]))

    %% Test Vars
    ylabel(ylab)
    xlabel(xlab)
    %% end - Test vars


    slab = {};
    for i = 1:numel(states)
        subplot(1,numel(states)+1,i+1),hold on
        ind = Trial.stc{states(i)};
        slab{end+1} = ind.label; 
        N = hist2([var1(ind),var2(ind)],v1edgs,v2edgs);
        imagesc(v1edgs,v2edgs,N',[0,mean(N(nniz(N(:))))*5])
        ylabel(ylab)
        xlabel(xlab)
        xlim(v1edgs([1,end]))
        ylim(v2edgs([1,end]))


    end
    suptitle(['JPDF ' Trial.filebase ])

    reportfig(fullfile(Trial.path.data,'figures'),figh,'jpdf_ang_height',[],Trial.filebase,100)

end
