% req20210503.m
%
%
% 
%

MjgER2016_load_data();

t = 20;

cluIds = units{t};

rmap = plot(pfs,25);


[mxr,mxp] = maxRate(pfs,cluIds);

% Multi-field placecells require an immutable space-theta-phase code to generate an unambiguous
% representation of space. 

% External input on top of topological constraints of the network 
%     - behavior state
%     - head direction???
%
% Visualization
%
% A trajectory extends outward 

% CA1 gaussian receptive field. 
% CA3 attractor coordinate.


fet = fet_dxy(Trials{t});

pfs = {};
for s = 1:numel(states);
state     = states{s};%'theta-groom-sit';
overwrite = false;

pargs = get_default_args('MjgER2016','MTAApfs','struct');
pargs.units = units{t};
pargs.states = state;
pargs.overwrite = overwrite;
pargs.boundaryLimits   = [-pi,pi;    ...
                          -500,500;  ...
                          -500,500];
pargs.binDims          = [0.2,20,20];
pargs.SmoothingWeights = [3,3,3];
pargs.numIter = 1;
pargs.halfsample = false;
pargs.xyzp = fet;
pargs.compute_pfs = @PlotPFCirc;

pfsArgs = struct2varargin(pargs);
pfs{s} = MTAApfs(Trials{t},pfsArgs{:});
end

alpha = linspace(-pi,pi,1000);
hfig = figure();
u = units{t}(2);
while u ~= -1 
    clf(hfig);
for s = 1:numel(states);    
    rmap = plot(pfs{s},u);
    mrmap = sq(mean(rmap,'omitnan'));
    mrmap(~nniz(mrmap(:))) = 0;

    mxi = LocalMinima2(-mrmap,1,200);
    mxr = [];
    if isempty(mxi); continue; end;
    for m = size(mxi,1),
        mxr(m) = mrmap(sub2ind([50,50],mxi(m,1),mxi(m,2)));
    end
    clim = [0,mxr];
    subplot2(5,7,1,s);    
        hold(hfig.CurrentAxes,'on');
        hp = pcolor(mrmap');
        hp.EdgeColor='none';        
        axis('xy');
        plot(mxi(1),mxi(2),'*m');
        axis('tight');
        title(states{s});
        caxis(clim);


    subplot2(5,7,2,s);
        hp = pcolor(rmap(:,:,mxi(2))');
        hp.EdgeColor='none';        
        axis('xy');
        title(num2str(u))
        caxis(clim);        
        
    subplot2(5,7,3,s);
        hp = pcolor(sq(rmap(:,mxi(1),:))');
        hp.EdgeColor='none';
        axis('xy');
        title(['max rate: ',num2str(mxr)]);        
        caxis(clim);        
        
        
        myaw = nan([50,50]);
        mr   = nan([50,50]);
        mrt  = nan([50,50]);
        for x = 1:50,
            for y = 1:50,
                cpx = rmap(:,x,y).*exp(-i*pfs{s}.adata.bins{1})./sum(rmap(:,x,y),'omitnan');
                if mrmap(x,y) >0.5,
                    myaw(x,y) = angle(sum(cpx,'omitnan'));
                    mr(x,y) = abs(sum(cpx,'omitnan'));                
                end
            end
        end


    subplot2(5,7,4,s);
        hold(hfig.CurrentAxes,'on');
        hp = pcolor(myaw');
        hp.EdgeColor='none';
        axis('xy'); 
        plot(mxi(1),mxi(2),'*m');
        axis('tight');
        scatter(24*cos(alpha)+25.5,24*sin(alpha)+25.5,10,alpha,'Filled');        
        colormap(hfig.CurrentAxes,'HSV');
        
        
        
    subplot2(5,7,5,s);
        hp = pcolor(mr');
        hp.EdgeColor='none';
        axis('xy'); 
        colormap(hfig.CurrentAxes,'JET');
        caxis(hfig.CurrentAxes,[0,1]);
      
end

    u = figure_controls(hfig,u,units{t});

end
%figure END



figure,imagesc(sq(rmap(:,33,:))');axis(xy);


figure,imagesc(sq(mean(rmap,'omitnan'))');axis('xy');