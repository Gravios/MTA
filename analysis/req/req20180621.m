

MjgER2016_load_data();


% LOAD all spk objects
% LOAD all lfp

%states = {'theta-groom-sit','rear&theta','hloc&theta','hpause&theta','lloc&theta','lpause&theta'};
states = {'theta','loc','pause','rear','hloc','hpause','lloc','lpause','groom','sit'};
distThresh = 300;

spkmap = [];
spkdrz = [];
spkpch = [];
spkhgt = [];
spksvd = [];
spkphz = [];
spkvxy = [];
nspkstc = [];
spkfbr = [];
snnpksvd = [];
spkavl = [];
spktrans = repmat({[]},[2,numel(states)]);

for tind = 1:numel(Trials);
    Trial = Trials{tind};
    unitSubset = units{tind};
    xyz = preproc_xyz(Trial,'trb');
    xyz.filter('ButFilter',3,2.5,'low');
    vxy = xyz.vel({'spine_lower','nose','hcom'},[1,2]);

    fet = fet_HB_pitchB(Trial);
    fsvd = fet_svd(Trial);
    fbr = fet_bref(Trial);

    xyz = preproc_xyz(Trial,'trb');            
    xyz.filter('ButFilter',3,1.5,'low');    
    ang = create(MTADang,Trial,xyz);

    avl = MTADfet.encapsulate(Trial,...
                              [circ_dist(circshift(ang(:,'head_back','head_front',1),-3),...
                                        circshift(ang(:,'head_back','head_front',1),3)).*xyz.sampleRate,...
                               circ_dist(circshift(ang(:,'spine_lower','spine_upper',1),-3),...
                                        circshift(ang(:,'spine_lower','spine_upper',1),3)).*xyz.sampleRate],...
                              ang.sampleRate,...
                              'angular velocity',...
                              'angvel',...
                              'v');
    
    % LOAD behavioral states
    stc = Trial.stc.copy();
    % LOAD units
    spk = Trial.load('spk',xyz.sampleRate,'gper',unitSubset,'deburst');
    pft = pfs_2d_theta(Trial,unitSubset);
    drz = compute_drz(Trial,unitSubset,pft);
    ddz = compute_ddz(Trial,unitSubset,pft);
    lfp = Trial.load('lfp',sessionList(tind).thetaRef);
    lfp.resample(xyz);
    phz = lfp.phase([5,13]);
    
    stcm = stc2mat(stc,xyz,states);
    
    ststrans = {};
    for state = 1:numel(states),
        samples = round(xyz.sampleRate/2);
        stsper = [Trial.stc{states{state}}];
        stsdiff = [[stsper.data(:,2);size(xyz,1)]-[0;stsper.data(:,1)]]';
        for s = 1:2,
            trans = stsper.data(:,s)';
            trans(stsdiff(s:end-2+s)<=samples) = [];
            trans((trans-samples)<=0|(trans+samples)>=size(xyz,1)) = [];
            ststrans{s,state} = xyz.copy('empty');
            ststrans{s,state}.data = nan([size(xyz,1),1]);
            
            for transition = trans,
                ststrans{s,state}.data(transition-samples:transition+samples) = ...
                    linspace(-1,1,2*samples+1);
            end
        end
    end


    for unit = 1:numel(unitSubset);

        % GET unit index    
        % GET spikes times of unit    
        % REMOVE spikes outside of position acquisition periods
        % REMOVE spikes outside of distance threshold
        res = spk(unitSubset(unit));
        res(res>size(drz,1))=[];
        res (abs(ddz(res, unit))>=distThresh) = [];            

        % GET direction rate zone (DRZ) values at times of spikes ( see Huxter(2008) )
        % GET phase values at times of spikes   
        % IGNORE spikes where drz or phase are nans
% $$$         for state = 1:numel(states);
% $$$             s = [stc{states{state},xyz.sampleRate}];
% $$$             sres = res(WithinRanges(res,s.data));
% $$$             spkdrz{state} = cat(1,spkdrz{state},drz(sres,unit));
% $$$             spkpch{state} = cat(1,spkpch{state},ang(sres,'spine_middle','spine_upper',2));
% $$$             spkhgt{state} = cat(1,spkhgt{state},xyz(res,'head_front',3));
% $$$             spkphz{state} = cat(1,spkphz{state},phz(sres,spk.map(unitSubset(unit)==spk.map(:,1),2)));
% $$$         end
        spkmap = cat(1,spkmap,[ones([numel(res),1]).*tind,...
                               ones([numel(res),1]).*unitSubset(unit)]);
        spkstc = cat(1,spkstc,stcm(res,:));
        spkpch = cat(1,spkpch,fet(res,:));        
        spkvxy = cat(1,spkvxy,vxy(res,:));        
        spkfbr = cat(1,spkfbr,fbr(res,:));        
        spksvd = cat(1,spksvd,fsvd(res,:));
        spkavl = cat(1,spkavl,avl(res,:));        
        spkhgt = cat(1,spkhgt,xyz(res,'nose',3));
        spkdrz = cat(1,spkdrz,drz(res,unit));
        spkphz = cat(1,spkphz,phz(res,spk.map(unitSubset(unit)==spk.map(:,1),2)));        
        for state = 1:numel(states),
            for s = 1:2,
                spktrans{s,state} = cat(1,spktrans{s,state},ststrans{s,state}(res));
            end
        end
    end
end

figure,plot3(spkdrz{end},spktrans{1,3},spkphz{end},'.');






% LOAD behavior scores
if ~exist('pfd','var'), [pfd,tags,eigVec,eigVar,eigScore,validDims,unitSubsets,unitIntersection,zrmMean,zrmStd] = req20180123_ver5(Trials);  end
numComp = size(eigVec{1},2);
pfindex = 1;
MjgER2016_load_bhv_erpPCA_scores();
% output:
%    fsrcz
%    FSrC
%    rmaps
sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);
cc = FSrC(:,[2,1,3])+0.75;
cc(~sigUnits,:) = repmat([0.75,0.75,0.75],[sum(~sigUnits),1]);
figure();
scatter3(FSrC(:,1),FSrC(:,2),FSrC(:,3),20,cc,'filled');
xlim([-4,4]);
ylim([-4,4]);
zlim([-4,4]);
% $$$ sp(end).YTickLabel = {};
% $$$ sp(end).XTickLabel = {};
box('on');
hold('on');
set(gca,'XColor',[0.75,0.75,0.75],'YColor',[0.75,0.75,0.75]);
cluMap = [20,74;...
          20,79;...
          20,104];
cluSessionSubset = cluSessionMap(unitSubsets{pfindex},:);
for u = cluMap'
    uind = find(ismember(cluSessionSubset,u','rows'));
    scatter3(FSrC(uind,1),FSrC(uind,2),FSrC(uind,3),40,[0,0,0]);    
end




figure();
binDims = [20,20];

binsTrans = linspace(-pi/2,pi/2,binDims(2));
binsTransInd = discretize(spkpch(:,2),binsTrans);
labelTrans = 'Body Pitch (rad)';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_bpitch';
 

binsTrans = linspace(-10,10,binDims(2));
binsTransInd = discretize(spksvd(:,2),binsTrans);
labelTrans = 'svdpc2';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_fsvd';


binsTrans = linspace(-pi/2,pi/2,binDims(2));
binsTransInd = discretize(spkpch(:,2),binsTrans);
labelTrans = 'Head Pitch (rad)';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_hpitch';

binsTrans = linspace(20,300,binDims(2));
binsTransInd = discretize(spkhgt,binsTrans);
labelTrans = 'Height (mm)';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_height';

binsTrans = linspace(-1,1.8,binDims(2));
binsTransInd = discretize(log10(spkvxy(:,2)),binsTrans);
labelTrans = 'Head Speed log10(cm/s)';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_hspeed';


binsTrans = linspace(-1,1.8,binDims(2));
binsTransInd = discretize(log10(spkvxy(:,1)),binsTrans);
labelTrans = 'Body Speed log10(cm/s)';
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = 'drz_x_bspeed';



sind = 9;
binsTrans = linspace(-1,1,binDims(2));
binsTransInd =  discretize(spktrans{1,sind},binsTrans);
labelTrans = [states{sind},' state onset'];
binsDrz = linspace(-1,1,binDims(1));
binsDrzInd  =  discretize(spkdrz,binsDrz);
labelDRZ = 'DRZ';
saveLabel = ['drz_x_',states{sind},'ON'];

              

%sigUnitsBhv = true([size(FSrC,1),1]);
%sigUnitsBhv = any(FSrC(:,[1,3])>=-0.5,2);
sigUnitsBhv = any(FSrC(:,[1,3])>=-0.5,2);
%sigUnitsBhv = any(FSrC(:,[2])<0,2);


ind = nniz(binsTransInd)&nniz(binsDrzInd)           ...
      & spkstc(:,1)                                     ... theta
      &~spkstc(:,9)                                     ... not groom
      &~spkstc(:,10)                                    ... not sit
      &~spkstc(:,4)                                     ... not rear
      & ismember(spkmap,cluSessionSubset(sigUnitsBhv,:),'rows');
A = accumarray([binsDrzInd(ind),binsTransInd(ind)],spkphz(ind),[numel(binsDrz)-1,numel(binsTrans)-1],@circ_mean);
S = accumarray([binsDrzInd(ind),binsTransInd(ind)],spkphz(ind),[numel(binsDrz)-1,numel(binsTrans)-1],@circ_std);
C = accumarray([binsDrzInd(ind),binsTransInd(ind)],ones([sum(ind),1]),[numel(binsDrz)-1,numel(binsTrans)-1],@sum);


clf();
subplot(131);imagesc(binsDrz,binsTrans,A');axis('xy');colormap(gca,'hsv');
xlabel(labelDRZ);ylabel(labelTrans);
cax = colorbar();ylabel(cax,'Mean Theta Phase');
subplot(132);imagesc(binsDrz,binsTrans,S');colorbar();axis('xy');colormap(gca,'default');
xlabel(labelDRZ);ylabel(labelTrans);
cax = colorbar();ylabel(cax,'STD Theta Phase');
subplot(133);imagesc(binsDrz,binsTrans,C');colorbar();axis('xy');colormap(gca,'default');
xlabel(labelDRZ);ylabel(labelTrans);
cax = colorbar();ylabel(cax,'Count');

hax = findobj(gcf,'Type','Axes');
af(@(h) set(h,'Units','centimeters'),  hax);
af(@(h) set(h,'Position',[h.Position(1:2),4,4]),  hax);

print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['pp_pop_',saveLabel,'.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['pp_pop_',saveLabel,'.png']]);






%sigUnits = any(abs(fsrcz(:,1:3))>=1.96,2);



figure();
for state = 1:numel(states)
    subplot(1,numel(states),state);
    hold('on');
    plot(spkdrz{state},spkphz{state},'.b');
    plot(spkdrz{state},spkphz{state}+2*pi,'.b');
end



%sigUnitsBhv = ismember(spkmap,cluSessionSubset(all(abs(fsrcz(:,1:3))<1.96,2),:),'rows');

%sigUnitsBhv = ismember(spkmap,cluSessionSubset(any(fsrcz(:,1)>1,2),:),'rows');
figure();
clf();
for state = 1:numel(states)
subplot(1,numel(states),state);
ind = logical(spkstc(:,1))&logical(spkstc(:,state));
    hold('on');
    hist2([[spkdrz(ind);spkdrz(ind)],...
           rad2deg([spkphz(ind);spkphz(ind)+2*pi])],...
          linspace(-1,1,50),rad2deg(linspace(-pi,3*pi,50)));
    axis('tight');
    colormap('jet');
    title(states{state});
    xlabel('drz');
    ylabel('phaze (rad)');    

    hax = gca();
    hax.YTick = [-180,-90,0,90,180,270,360,450,540];
    hax.Units = 'centimeters';
    hax.Position = [hax.Position(1:2),3,3];    
end
print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['pp_pop_state.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['pp_pop_state.png']]);



%sigUnitsBhv = ismember(spkmap,cluSessionSubset(all(abs(fsrcz(:,1:3))<1.96,2),:),'rows');
sigUnitsBhv = ismember(spkmap,cluSessionSubset(any(FSrC(:,[2])>0,2),:),'rows');
sigUnitsBhv = ismember(spkmap,cluSessionSubset(any(FSrC(:,[2])>0,2),:),'rows');
%sigUnitsBhv = ismember(spkmap,cluSessionSubset(any(fsrcz(:,1)>1,2),:),'rows');
figure();
clf();
ind = sigUnitsBhv&logical(spkstc(:,1))&logical(spkstc(:,4));
%ind = sigUnitsBhv;%&logical(spkstc(:,1));
%for state = 1:numel(states)
%subplot(1,numel(states),state);
    hold('on');
    hist2([[spkdrz(ind);spkdrz(ind)],...
           rad2deg([spkphz(ind);spkphz(ind)+2*pi])],...
          linspace(-1,1,50),rad2deg(linspace(-pi,3*pi,50)));
% $$$     hist2([[spkpch(ind,2);spkpch(ind,2)],...
% $$$            rad2deg([spkphz(ind);spkphz(ind)+2*pi])],...
% $$$           linspace(-0.6,1.5,50),rad2deg(linspace(-pi,3*pi,50)));
% $$$     hist2([[spkhgt(ind,1);spkhgt(ind,1)],...
% $$$            rad2deg([spkphz(ind);spkphz(ind)+2*pi])],...
% $$$           linspace(10,250,50),rad2deg(linspace(-pi,3*pi,50)));
% $$$     hist2([[spkhgt(ind,1);spkhgt(ind,1)],...
% $$$            rad2deg([spkphz(ind);spkphz(ind)+2*pi])],...
% $$$           linspace(10,250,50),rad2deg(linspace(pi/2,5/2*pi,50)));
% $$$     hist2([[spktrans{1,4}(ind);spktrans{1,4}(ind)],...
% $$$            rad2deg([spkphz(ind);spkphz(ind)+2*pi])],...
% $$$           linspace(-1,1,30),rad2deg(linspace(-pi,3*pi,30)));
    
    axis('tight');
    colormap('jet');
    title(states{state});
    xlabel('drz');
    ylabel('phaze (rad)');    

    hax = gca();
    hax.YTick = [0,90,180,270,360];
    hax.Units = 'centimeters';
    hax.Position = [hax.Position(1:2),3,3];
    %end






Trial = Trials{20};
unitSubset = [20,21,25,35,79,83,119];        
xyz = Trial.load('xyz');
% LOAD behavioral states
stc = Trial.stc.copy();
% LOAD units
spk = Trial.load('spk',xyz.sampleRate,'theta-groom-sit',unitSubset,'deburst');
pft = pfs_2d_theta(Trial,unitSubset);
drz = compute_drz(Trial,unitSubset,pft);
ddz = compute_ddz(Trial,unitSubset,pft);
lfp = Trial.load('lfp',sessionList(tind).thetaRef);
lfp.resample(xyz);
phz = lfp.phase([5,13]);

figure();

for unit = 1:numel(unitSubset);
    % GET unit index    
    % GET spikes times of unit    
    % REMOVE spikes outside of position acquisition periods
    % REMOVE spikes outside of distance threshold
    res = spk(unitSubset(unit));
    res(res>size(drz,1))=[];
    res (abs(ddz(res, unit))>=distThresh) = [];            

    % SKIP fitting parameters if too few spikes
    if numel(res)<10,  continue;  end    
    
    % GET direction rate zone (DRZ) values at times of spikes ( see Huxter(2008) )
    % GET phase values at times of spikes   
    % IGNORE spikes where drz or phase are nans
    clf();
    subplot(1,numel(states)+1,1);        
    pft.plot(unitSubset(unit),'mean',false,[],true,'interpPar',interpParPfs,'colorMap',@jet);
    hax = gca();
    hax.Units = 'centimeters';
    hax.Position = [hax.Position(1:2),3,3];
    
    for state = 1:numel(states);
        s = [stc{states{state},xyz.sampleRate}];
        sres = res(WithinRanges(res,s.data));
        spkdrz = drz(sres,unit);
        spkphz = phz(sres,spk.map(unitSubset(unit)==spk.map(:,1),2));
        spkgid = ~
        isnan(spkdrz)&~isnan(spkphz);
        subplot(1,numel(states)+1,state+1);
        plot([spkdrz(spkgid);spkdrz(spkgid)],...
             rad2deg([spkphz(spkgid);spkphz(spkgid)+2*pi]),...
             '.','MarkerSize',2);
        hax = gca();
        hax.YTick = [-180,-90,0,90,180,270,360];
        hax.Units = 'centimeters';
        hax.Position = [hax.Position(1:2),3,3];
        axis('tight');

    end

    print(gcf,'-depsc2',...
          ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
           ['pp_example_',Trial.filebase,'_unit_',unitSubset(unit),'.eps']]);
    print(gcf,'-dpng',...
          ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
           ['pp_example_',Trial.filebase,'_unit_',unitSubset(unit),'.png']]);


end








trialIndex = 20;

Trial = Trials{trialIndex}; 
unitSubset = units{trialIndex};


% COMPUTE Ndim placefield
xyzp = []
pfsArgs = struct('units',              unitSubset,                           ...
                 'states',             'theta-groom-sit',                    ...
                 'overwrite',          false,                                ...
                 'tag',                '',                                   ...
                 'binDims',            [ 100, 100,0.4,0.4],                  ...
                 'SmoothingWeights',   [0.8,0.8,0.8,0.8],                    ...
                 'type',               'xyhb',                               ...
                 'spkShuffle',         false,                                ...
                 'posShuffle',         false,                                ...
                 'numIter',            1000,                                 ...
                 'xyzp',               xyzp,                                 ...
                 'boundaryLimits',     [-500,500;-500,500;-2,2;-2,2],        ...
                 'bootstrap',          false,                                ...
                 'halfsample',         true                                  ...
                 );
pfsArgs = struct2varargin(pfsArgs);
pfs = MTAApfs(Trial,pfsArgs{:});

stc = Trial.load('stc','msnn_ppsvd_raux');
xyz = resample(preproc_xyz(Trial,'trb'),10);
fet = fet_HB_pitchB(Trial,10);
ufr = Trial.ufr.copy;
ufr = ufr.create(Trial,xyz,'gper',unitSubset,0.75,true);
