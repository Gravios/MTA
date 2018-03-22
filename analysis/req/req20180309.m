
% req20180309 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Distance vs drz restricted behavior field egenvector scores
%               
%  Bugs: NA


figDir = '/storage/gravio/figures/analysis/placefields_nonSpatialFeatures';



sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);


Trials = af(@(s) MTATrial.validate(s), sessionList);

[pfd,tags,eigVec,eigVar,eigScore,unitSubsets,unitIntersection] = req20180123_ver5(Trials);


units = cf(@(T)  select_placefields(T),  Trials); 
units = req20180123_remove_bad_units(units);

cluMap = [];
for u = 1:numel(units)
    cluMap = cat(1,cluMap,[u*ones([numel(units{u}),1]),units{u}(:)]);
end

pft = cf(@(T,u)  pfs_2d_theta(T,u,'overwrite',false),  Trials,units);

[mxr,mxp] = cf(@(p,u) p.maxRate(u),pft,units);

uMaxPos = cat(1,mxp{:});

% PLOT radial distance from maze center againts "rear" score

pfStates = {}
pfStates{1} = {'low','rear','high'};
pfStates{2} = {'low','high'};
pfStates{3} = {'low','high'};

pfindex = 1;
for s = 1:numel(pfStates{pfindex}),
    hfig = figure(382939); clf();
    hfig.Units = 'centimeters';
    hfig.Position = [0.5,0.5,10,6];
    hfig.PaperPositionMode = 'auto';
    hax = subplot(1,1,1);
    plot(sqrt(sum(uMaxPos(unitSubsets{pfindex}(ismember(unitSubsets{pfindex},unitIntersection)),:).^2,2)),...
         eigScore{pfindex}(ismember(unitSubsets{pfindex},unitIntersection),s),...                           
         '.','MarkerSize',2);
    xlabel('radial distance (mm)');
    ylabel([pfStates{pfindex}{s},' behavior score']);
    title({'radial distance from maze center',...
           'vs ',pfStates{pfindex}{s},' behavior score'});
    
    xlim([0,450]);
    ylim([-2,6]);    

    af(@(h) set(h,'Units','centimeters'),            hax);    
    af(@(h) set(h,'Position',[4,1,1.5,1.5]), hax);

    figName = ['erpPCA_',tags{pfindex},'_mazeDist_',pfStates{pfindex}{s},'Score'];
    print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
    print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));
end


hfig = figure(382940); 
hfig.Units = 'centimeters';
hfig.Position = [0.5,0.5,10,6];
hfig.PaperPositionMode = 'auto';
hold('on')
clf();
pfindex = 1;
scl ='brg';
hax = subplot(1,1,1);
for s = 1:numel(pfStates{pfindex}),
mRadialDist =sqrt(sum(uMaxPos(unitSubsets{pfindex}(ismember(unitSubsets{pfindex},unitIntersection)),:).^2,2));
subEigScore = eigScore{pfindex}(ismember(unitSubsets{pfindex},unitIntersection),s);
[radialEcdf,radDist,recdfUp,recdfDown] = ecdf(mRadialDist(subEigScore>1));
boundedline(radDist,radialEcdf,[radialEcdf-recdfDown,recdfUp-radialEcdf],['-',scl(s)],'alpha')
end
title(['erpPCA scores ',strjoin(pfStates{pfindex},'_'),' ',tags{pfindex}])
af(@(h) set(h,'Units','centimeters'),    hax);    
af(@(h) set(h,'Position',[4,1,1.5,1.5]), hax);       
figName = ['erpPCA_',tags{pfindex},'_mazeDist_',strjoin(pfStates{pfindex},'-'),'Score_ecdf'];
print(hfig,'-depsc2',fullfile(figDir,[figName,'.eps']));        
print(hfig,'-dpng',  fullfile(figDir,[figName,'.png']));






% PLOT radial distance from maze center against "low loc" score
figure,
plot(sqrt(sum(uMaxPos(unitSubsets{1}(ismember(unitSubsets{1},unitIntersection)),:).^2,2)),...   
     eigScore{1}(ismember(unitSubsets{1},unitIntersection),1),...
     '.');



figure,
histc(sqrt(sum(uMaxPos(unitSubsets{1}(ismember(unitSubsets{1},unitIntersection)),:).^2,2)),...   % radial distance from maze center
     eigScore{1}(ismember(unitSubsets{1},unitIntersection),2),...                           % "rear" score         
     '.');



figure,
scatter(eigScore{1}(ismember(unitSubsets{1},unitIntersection),1),...
        eigScore{1}(ismember(unitSubsets{1},unitIntersection),3),...
        10,...
        sqrt(sum(uMaxPos(unitSubsets{1}(ismember(unitSubsets{1},unitIntersection)),:).^2,2)),...
        'filled');



