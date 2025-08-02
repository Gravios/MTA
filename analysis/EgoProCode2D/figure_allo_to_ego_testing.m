

configure_default_args();
MjgER2016_load_data();
%EgoProCode2D_load_data();


% COMPUTE the size of primary egofields

%%%<<< Load general variables
sampleRate = 250;
%headCenterCorrection = [-25,-8];
pfsState = 'theta-groom-sit-rear';
hbaBinEdges = -1.5:0.6:1.5;
xyz = cf(@(t) preproc_xyz(t,'trb'),             Trials);
      cf(@(x) x.filter('ButFilter',3,30,'low'), xyz);    
      cf(@(x) x.resample(sampleRate),           xyz);
spk = cf(@(t,u) t.load('spk',sampleRate,'gper',u,'deburst'),Trials,units);    
pft = cf(@(t,u)  pfs_2d_theta(t,u),  Trials, units);
%%%>>>


tind = 18;
overwrite = false;
%%%<<< (pfe) ego ratemap
pfe = cf(@(t,u,x,s,p)                                          ... Egocentric ratemap .
         compute_ego_ratemap(t,u,x,s,p,'overwrite',overwrite), ...
             Trials(tind),                                           ... MTATrial
             units(tind),                                            ... Unit subset, placefields away from the maze walls
             xyz(tind),                                              ... MTADxyz object, head position
             spk(tind),                                              ... MTASpk object, spike time and id collection 
             pft(tind)                                               ... MTAApfs object, theta state placefields 
);
%%%>>>

%%%<<< (pfet) egothp ratemap
pfet = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase.
         compute_egothp_ratemap(t,u,x,s,p,'overwrite',overwrite),...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind),                           ... MTASpk object, spike time and id collection 
             pft(tind)                            ... MTAApfs object, theta state placefields 
);
%%%>>>

%%%<<< (pfs) egohba ratemap
pfs = cf(@(t,u,x,s,p)                       ... Egocentric ratemap | theta phase , head body angle.
         compute_egohba_ratemap(t,u,x,s,p,'overwrite',overwrite),   ...
             Trials(tind),                        ... MTATrial
             units(tind),                         ... Unit subset, placefields away from the maze walls
             xyz(tind),                           ... MTADxyz object, head position
             spk(tind),                           ... MTASpk object, spike time and id collection 
             pft(tind)                            ... MTAApfs object, theta state placefields 
);
% $$$ pfsh = cf(@(t,u,x,s,p)                      ... Egocentric ratemap | theta phase , shuffled head body angle.
% $$$          compute_egohba_ratemap_shuffled(t,u,x,s,p,'overwrite',overwrite),   ...
% $$$              Trials,                        ... MTATrial
% $$$              units,                         ... Unit subset, placefields away from the maze walls
% $$$              xyz,                           ... MTADxyz object, head position
% $$$              spk,                           ... MTASpk object, spike time and id collection 
% $$$              pft                            ... MTAApfs object, theta state placefields 
% $$$ );
%%%>>>


% $$$ 
% $$$ 
% $$$ figure();
% $$$ for u = 1:numel(units{tind});
% $$$     subplot2(5,2,1,1);
% $$$         plot(pft{tind},units{tind}(u),1,'colorbar',[],true);
% $$$         title(num2str(units{tind}(u)));
% $$$     subplot2(5,2,2,1);
% $$$         plot(pfe{1},units{tind}(u),1,'colorbar',[],false);
% $$$         mrate = max(cell2mat(cf(@(p) max(p.data.rateMap(:,p.data.clu==units{t}(u))), pfet{tind})));
% $$$     for p = 1:5
% $$$         subplot2(5,2,p,2);
% $$$             plot(pfet{1}{p},units{tind}(u),1,'colorbar',[0,mrate],false,'colorMap',@jet);
% $$$     end
% $$$     waitforbuttonpress();
% $$$ end




figure();
for u = 1:numel(units{tind});
    subplot2(5,7,1,1);
        plot(pft{tind},units{tind}(u),1,'colorbar',[],true);
        title(num2str(units{tind}(u)));
    subplot2(5,7,2,1);
        plot(pfe{1},units{tind}(u),1,'colorbar',[],false,'flipAxesFlag',true);
        mrate = max(cell2mat(cf(@(p) max(p.data.rateMap(:,p.data.clu==units{tind}(u))), pfet{1})));
        xlim([-300,300]);ylim([-300,300]);
    for p = 1:5
        subplot2(5,7,p,2);
            plot(pfet{1}{6-p},units{tind}(u),1,'colorbar',[0,mrate],false,'flipAxesFlag',true,'colorMap',@jet);
            Lines(0,[],'m');
            Lines([],0,'m');
            xlim([-300,300]);ylim([-300,300]);            
    end
    for p = 1:5
        for h = 1:5
        subplot2(5,7,p,h+2);
            plot(pfs{1}{6-p,h},units{tind}(u),1,'colorbar',[0,mrate],false,'flipAxesFlag',true,'colorMap',@jet);
            Lines(0,[],'m');
            Lines([],0,'m');
            xlim([-300,300]);ylim([-300,300]);
        end
    end
    waitforbuttonpress();
end

    
    
%[hfig,fig,fax,sax] = set_figure_layout(figure(666001),'A4','portrait',[],1.5,1.5,0.05,0.05);
