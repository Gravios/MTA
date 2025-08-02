MjgER2016_load_data();

if ~exist('pfd','var'), [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials);  end
numComp = size(eigVec{1},2);
pfindex = 1;
MjgER2016_load_bhv_erpPCA_scores();
% output:
%    fsrcz
%    FSrC
%    rmaps
clear('pfd','tags','eigVec','eigVar','eigScore','validDims','zrmMean','zrmStd',...
      'clu','tlu','rind','D','LR','FSCFr','rsMean','rsStd','pfdShuffled','rmapsShuffledMean',...
      'rmapsShuffled','FSrM','FSrS','fsrsMean','fsrsStd','fsrsMean','fsrsStd');


hfig = figure();
% STATE - locomotion
nPart = numel(modelTraj.partitions)-1;

clf();
hfig.Units = 'centimeters';
hfig.Position = [0,0,30,5*nPart];
for v = 1:nPart
    y = nPart-v+1;
    % select data with in speed partition v
    ind = cind  &  sind  ...
                &  vxy(:,speedInd)>modelTraj.partitions(v)   ...
                &  vxy(:,speedInd)<modelTraj.partitions(v+1);
% REFERENCE trajectory coordinate system
    %decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
    decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),tvec(ind,:,:),2,[2,3]);
    
    % Anterior-Posterior axis
    subplot2(nPart,6,y,1);
    % JPDF phase X error
    hist2([[decError(:,1);decError(:,1)],...
           [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-300,300,40),...
          linspace(-pi,pi*3,30)); 
    line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300])-2*pi,'Color','m');    
    line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300]),     'Color','m');
    line([-300,300],polyval(modelTraj.parameters(v,:),[-300,300])+2*pi,'Color','m');
    title(['rho: ',num2str(modelTraj.rho(v))]);
    xlabel('mm');
    ylabel('theta phase');
    Lines(0,[],'k');
% medial-lateral axis
    subplot2(nPart,6,y,2);
    hist2([[decError(:,2);decError(:,2)],...
           [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-300,300,40),...
          linspace(-pi,pi*3,30)); 
    xlabel('mm');
    ylabel('theta phase');
    Lines(0,[],'k');
% direction
    subplot2(nPart,6,y,3);
    hist2([[atan2([decError(:,2);decError(:,2)],...
                  [decError(:,1);decError(:,1)])],...
           [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-pi,pi,40),...
          linspace(-pi,pi*3,30)); 
    xlabel('yaw');
    ylabel('theta phase');
    Lines(0,[],'k');
% REFERENCE head coordinate system
    %decError = multiprod(posEstCom(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);    
    decError = multiprod(decPos(ind,[1,2])-sq(xyz(ind,'hcom',[1,2])),hvec(ind,:,:),2,[2,3]);    
    subplot2(nPart,6,y,4);
    hist2([[decError(:,1);decError(:,1)],...
          [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-300,300,40),...
          linspace(-pi,pi*3,30)); 
    line([-300,300],polyval(modelHead.parameters(v,:),[-300,300])-2*pi,'Color','m');
    line([-300,300],polyval(modelHead.parameters(v,:),[-300,300]),'Color','m');
    line([-300,300],polyval(modelHead.parameters(v,:),[-300,300])+2*pi,'Color','m');
    
    title(['rho: ',num2str(modelHead.rho(v))]);
    xlabel('mm');
    ylabel('theta phase');
    Lines(0,[],'k');
% lateral
    subplot2(nPart,6,y,5);
    hist2([[decError(:,2);decError(:,2)],...
           [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-300,300,40),...
          linspace(-pi,pi*3,30)); 
    xlabel('mm');
    ylabel('theta phase');
    Lines(0,[],'k');
% yaw
    subplot2(nPart,6,y,6);
    hist2([atan2([decError(:,2);decError(:,2)],...
                 [decError(:,1);decError(:,1)]),...
          [phz(ind,1);phz(ind,1)+pi*2]],...
          linspace(-pi,pi,40),...
          linspace(-pi,pi*3,30));
    Lines(0,[],'k');
    xlabel('yaw');
    ylabel('theta phase');
end
