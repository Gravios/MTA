% req20200630 - Theta phase restricted egocentric placefield ratemaps
%    Tags: phase egocentric ratemap split by head body angle 
%    Status: active
%    Type: Analysis
%    Author: Justin Graboski
%    Project: MjgER2016
%    Description: Egocentric ratemaps decomposed by theta phase
%    Figure: none
%
%    processed:
%        compute_ego_ratemap
%        compute_egohba_ratemap
%        compute_egohba_ratemap_shuffled
%
%    ANALYSIS : 
%        Which behavioral variables affect phase precession?
%            HBA : head-body-angle
%            Which characteristics of phase precession change with respect to HBA?
%                The Forward and lateral 
%        compute_egohba_ratemap is segmented by theta phase
%        KNOWN : Lateralization of the egocentric rate map with increasing head-body angles on the ascending
%                phase of theta.
%        TEST  : shuffle the head-body-angles and compare the spatial information and rate weighted mean
%                position of each head-body-angle deflection bin. 
%        QUESTION : Is the spatial information significantly different in the shuffled condition?
%
%        DEFINITIONS : 
%            ERM : egocentric rate map
%            HBA : head-body-angle
%            ESI : egocentric spatial information
%            RWMP : rate-weigted-mean-position
%
%        VARIABLES : 
%            ESI of ERM-HBA
%            ESI of ERM-HBA(SHUFFLED)
%            RWMP of ERM-HBA
%            RWMP of ERM-HBA(SHUFFLED)
%            ZSCORE_ESI(ESI,ESI(shuffled))
%            ZSCORE_RWMP(RWMP,RWMP(shuffled))
%
%        PROCESS : ESI <- information <- rate maps <- MTAApfs objects
%
%        PLOT :  Joint Distribution of ZSCORE_ESI and ZSCORE_RWMP
%            SUBPLOT : for each HBA
%            

MjgER2016_load_data();

%%%<<< Define Unit Set
% SELECTED UNITS
units = {};
units{1}  = [ ]; %15, 96, 99,150,207                                % er01-20110719  CA3
units{2}  = [ 23, 86, 99,152,157,159]; %39,74, 75,                  % er01-20110721  CA3

units{3}  = [158,173];%95,                                          % ER06-20130612  CA1
units{4}  = [ 29, 31,107];%126,175,197                              % ER06-20130613  CA1
units{5}  = [ 31, 35, 67, 99, 121];%40, 48, 81,119,                 % ER06-20130614  CA1

units{6}  = [  4,  7, 10, 18, 25, 35, 37, 38, 49, 68,104];%60,      % Ed10-20140816  CA3
units{7}  = [  1, 10, 24, 33, 38, 57, 63, 64, 73, 82,105,108];      % Ed10-20140817  CA3

units{8}  = [ ];%                                                   % jg04-20120128  CA1
units{9}  = [ ];%23, 28];                                           % jg04-20120129  CA1
units{10} = [ ];% 1, 28];                                           % jg04-20120130  CA1
units{11} = [ ];%14, 15];                                           % jg04-20120131  CA1
units{12} = [ ];% 2];                                               % jg04-20120201  CA1

units{13} = [   ];                                                  % jg04-20120210  CA3
units{14} = [ 20];                                                  % jg04-20120211  CA3
units{15} = [   ]; %18                                              % jg04-20120212  CA3
units{16} = [  5];                                                  % jg04-20120213  CA3
units{17} = [   ]; %16, 20, 30, 32, 34, 44, 58, 61, 81 ,107];%  7,        % jg05-20120309  CA1
units{18} = [  9, 11, 21, 33, 35, 42, 49, 50, 54, 60, 63,...% jg05-20120310  CA1
              74, 75];                        % 18, 29, 78, 80
units{19} = [ 10, 13, 27, 28, 33, 50, 53, 63, 66, 67, 97,108,115,...% jg05-20120311  CA1 
             142,152,153,160,172];                % 73,94,139,              
units{20} = [ 20, 21, 25, 31, 35, 61, 72, 79, 80, 81, 85,...% jg05-20120312  CA1
             103,104,110,111,116,139,151];    %44, 52,109,138,
units{21} = [  6, 22, 24, 25, 37, 43, 61, 63, 77];% 33, 44, 62, 68,  % jg05-20120315  CA1

units{22} = [ 13, 18, 30, 41, 48, 56, 58, 61, 65];% 16,19(okwp),42, % jg05-20120316  CA1
units{23} = [ 10, 29, 40, 50, 54, 63, 72];        %35,              % jg05-20120317  CA1

% -- New Units --% 
units{24} = [ 18, 26, 48]; %  24,                                   % jg05-20120323  CA1  rear center 28, 52, 53
units{25} = [ 10, 29, 55]; %  12, 19,                               % jg05-20120324  CA1

units{26} = [  4, 43, 47, 70, 82,140,145,175,230];%87,155,167,211,  % ER06-20130624  CA3

units{27} = [ 39, 51, 83, 92, 98,100,102];                          % Ed10-20140815  CA3

units{28} = [ 91]; % 135,145,165,183,187                            % er01-20110722  CA3

% Define Unit Set
%%%>>>

%%%<<< Define inter-subject corrections and maps
rot = [0,0,                                         ... er01
       0,0,0,                                       ... ER06
       0,0,                                         ... Ed10
       0,0,0,0,0,0,0,0,0,                           ... jg04
       0.17, 0.17, 0.17, 0.17,0.17, 0.17, 0.17,     ... jg05
       0.17,0.17,0,0,0];% new units - jg05, jg05, ER06, Ed10, er01

brot =[0,0,                                         ... er01
       0,0,0,                                       ... ER06
       0,0,                                         ... Ed10
       0,0,0,0,0,0,0,0,0,                           ... jg04
       0, 0, 0, 0, 0, 0, 0,     ... jg05
       0, 0, 0, 0, 0];% new units - jg05, jg05, ER06, Ed10, er01

hbaCorrection = [-0.25,-0.25,                               ... er01
                 0.2,0.2,0.2,                       ... ER06
                 0,0,                               ... Ed10
                 0,0,0,0,0,0,0,0,0,                 ... jg04
                 -0.25.*ones([1,9]),                ... jg05
                0.2,0,-0.25]; % new units - jg05, jg05, ER06, Ed10, er01
hrlCorrection =[0,0,                                ... er01
                -0.27.*ones([1,3]),                 ... ER06
                -0.42,-0.42,                        ... Ed10
                -0.05.*ones([1,9]),                 ... jg04
                -0.48.*ones([1,9]),                 ... jg05
                -0.27,-0.42,0]; %                   ... ER06, Ed10, er01

phzCorrection = [pi*1.25,pi*1.25,                             ... er01
                 pi/2,pi/2,pi/2,                    ... ER06
                 pi/1.25,pi/1.25,                             ... Ed10
                 0,0,0,0,0,0,0,0,0,                 ... jg04                 
                 pi/4,pi/4,pi/4,pi/4,pi/4,pi/4,pi/4 ... jg05
                 pi/4,pi/4,pi,pi/1.25,pi*1.25]; % new units - jg05, jg05, ER06, Ed10, er01

anatomicalLocation = [0,0,                          ... er01
                      1,1,1,                        ... ER06
                      0,0,                          ... Ed10
                      0,0,0,0,0,1,1,1,1             ... jg04                      
                      1,1,1,1,1,1,1,                ... jg05
                      1,1,0,0,0];% new units - jg05, jg05, ER06, Ed10, er01                    
%%%>>>

%%%<<< Load general variables
sampleRate = 250;
headCenterCorrection = [-25,-8];
pfsState = 'theta-groom-sit-rear';
spkMode = 'deburst';
binPhzs = linspace(-pi,pi,6);
binPhzc = (binPhzs(1:end-1)+binPhzs(2:end))./2;

hbaBinEdges = -1.5:0.6:1.5;

thetaChan = [sessionList.thetaRefGeneral];

xyz = cf(@(t) preproc_xyz(t,'trb'),             Trials);
      cf(@(x) x.filter('ButFilter',3,30,'low'), xyz);    
      cf(@(x) x.resample(sampleRate),           xyz);

spk = cf(@(t,u) t.load('spk',sampleRate,'gper',u,'deburst'),Trials,units);    

pft = cf(@(t,u)  pfs_2d_theta(t,u,'pfsArgsOverride',...
                              struct('halfsample',false,'numIter',1)),           ...
                              Trials, units);
%%%>>>


%%%<<< ego ratemap
pfe = cf(@(t,u,x,s,p,rt,hc)                                    ... Egocentric ratemap given theta phase and head body angle.
         compute_ego_ratemap(t,u,x,s,p,rt,hc,headCenterCorrection,true),           ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection)                           ... head body angle correction (horizontal plane)
);
%%%>>>

%%%<<< egothp ratemap
pfet = cf(@(t,u,x,s,p,rt,hc,ch,pc)                             ... Egocentric ratemap given theta phase and head body angle.
         compute_egothp_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);
%%%>>>

%%%<<< egohba ratemap
pfs = cf(@(t,u,x,s,p,rt,hc,ch,pc)                              ... Egocentric ratemap given theta phase and head body angle.
         compute_egohba_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);

pfsh = cf(@(t,u,x,s,p,rt,hc,ch,pc)                             ... Egocentric ratemap given theta phase and head body angle.
         compute_egohba_ratemap_shuffled(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);
%%%>>>

%%%<<< ego hba hvl ratemap
pfl = cf(@(t,u,x,s,p,rt,hc,ch,pc)                              ... Egocentric ratemap given (TP, HBA, HVL gt 0)
         compute_egohbahvl_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed 
             num2cell(phzCorrection)                           ... theta phase offset
);

pflh = cf(@(t,u,x,s,p,rt,hc,ch,pc)                             ... Egocentric ratemap given (TP, HBA, HVL gt 0)
         compute_egohbahvl_ratemap_shuffled(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed 
             num2cell(phzCorrection)                           ... theta phase offset
);
%%%>>>

%%%<<< egohrl

pfr = cf(@(t,u,x,s,p,rt,hc,ch,pc)                              ... Egocentric ratemap given theta phase and head body angle.
         compute_egohrl_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hrlCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);

pfrh = cf(@(t,u,x,s,p,rt,hc,ch,pc)                             ... Egocentric ratemap given theta phase and head body angle.
         compute_egohrl_ratemap_shuffled(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hrlCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);
% egohrl
%%%>>>

pfs = cf(@(t,u,x,s,p,rt,hc,ch,pc)                              ... Egocentric ratemap given theta phase and head body angle.
         compute_egohba_ratemap(t,u,x,s,p,rt,hc,ch,pc,headCenterCorrection,true),   ...
             Trials,                                           ... MTATrial
             units,                                            ... Unit subset, placefields away from the maze walls
             xyz,                                              ... MTADxyz object, head position
             spk,                                              ... MTASpk object, spike time and id collection 
             pft,                                              ... MTAApfs object, theta state placefields 
             num2cell(rot),                                    ... head angle correction (horizontal plane)
             num2cell(hbaCorrection),                          ... head body angle correction (horizontal plane)
             num2cell(thetaChan),                              ... lfp channel from which phase is computed
             num2cell(phzCorrection)                           ... theta phase offset
);


%%%<<< ego tan rad ratemap
% SELECT only tangential trajectories
pfTan = cf(@(t,u,x,s,p,rt,hc,ch,pc)                            ... Egocentric ratemap given (TP, HBA, HVL gt 0)
           compute_egohbahvl_tanTraj_ratemap(t,u,x,s,p,rt,hc,ch,pc,true),...
               Trials,                                         ... MTATrial
               units,                                          ... Unit subset, placefields away from the maze walls
               xyz,                                            ... MTADxyz object, head position
               spk,                                            ... MTASpk object, spike time and id collection 
               pft,                                            ... MTAApfs object, theta state placefields 
           num2cell(rot),                                      ... head angle correction (horizontal plane)
           num2cell(hbaCorrection),                            ... head body angle correction (horizontal plane)
           num2cell(thetaChan),                                ... lfp channel from which phase is computed 
           num2cell(phzCorrection)                             ... theta phase offset
);

% SELECT only radial trajectories
pfRad = cf(@(t,u,x,s,p,rt,hc,ch,pc)                            ... Egocentric ratemap given (TP, HBA, HVL gt 0)
           compute_egohbahvl_radTraj_ratemap(t,u,x,s,p,rt,hc,ch,pc,true),...
               Trials,                                         ... MTATrial
               units,                                          ... Unit subset, placefields away from the maze walls
               xyz,                                            ... MTADxyz object, head position
               spk,                                            ... MTASpk object, spike time and id collection 
               pft,                                            ... MTAApfs object, theta state placefields 
           num2cell(rot),                                      ... head angle correction (horizontal plane)
           num2cell(hbaCorrection),                            ... head body angle correction (horizontal plane)
           num2cell(thetaChan),            ... lfp channel from which phase is computed 
           num2cell(phzCorrection)                             ... theta phase offset
);
%%%>>>

pfs = cf(@(t,u,x,s,p,r,rt,hc) ...
         compute_egohba_ratemap(t,u,x,s,p,rt,hc,ch),...
         Trials(18:23),units(18:23),xyz(18:23),spk(18:23),pft(18:23),...
         num2cell(rot(18:23)),...
         num2cell(hbaCorrection(18:23)),...
         num2cell([sessionList(18:23).thetaRefGeneral]) ...
);




figure,
for a = 1:5
    subplot2(5,2,a,1);
        plot(pfRad{20}{4},21,{':',':',a},'text',[0,10],false,[],true);
        xlim([-250,250]);
        ylim([-250,250]);
    subplot2(5,2,a,2);
        plot(pfTan{20}{4},21,{':',':',a},'text',[0,10],false,[],true);
        xlim([-250,250]);
        ylim([-250,250]);
end



figure,plot(pfe{20},25,1,'colorbar',[],false);


% GENERATE computational mask 
lims = {[-250,250],[-250,250]};
maskBinInd = cf(@(p,l) l(1) < p & p < l(2),  pfl{t}{1}.adata.bins([1:2]),lims);
mask = maskBinInd{1};
for dim = 2:numel(maskBinInd),
     mask = bsxfun(@and,mask,permute(maskBinInd{dim},[2:dim,1]));
end
maskDimVal = cf(@(p,i) p(i), pfl{t}{1}.adata.bins(1:2),maskBinInd);
maskBinPos = cell(size(maskDimVal));
[maskBinPos{:}] = ndgrid(maskDimVal{:});
maskBinPos = reshape(cat(ndims(maskBinPos)+1,maskBinPos{:}),[],ndims(maskBinPos));

% COMPUTE the rate weighted field postion
rweightedPos = {};
tid = {};
for t = 1:numel(Trials); 
    rweightedPos{t} = nan([numel(units{t}),...
                        2,...
                        numel(pfl{1}),...
                        pfl{1}{1}.adata.binSizes(end)]);
    for u = 1:numel(units{t}),
        for p = 1:numel(pfl{t}),
            rmap = plot(pfs{t}{p},units{t}(u));
            % Ratemap weighted position of egocetric placefield
            %     maskBinPos - coordinates of each bin
            %     ratemap normalized to 1
            al = 1:5;
            for a = 1:numel(al)
                rweightedPos{t}(u,:,p,a) = ...
                    sum(bsxfun(@times, ...
                               maskBinPos, ...                           
                               reshape(rmap(maskBinInd{:},a),[],1) ...
                               ./sum(reshape(rmap(maskBinInd{:},a),[],1),'omitnan')) ...
                        ,'omitnan');
            end
        end
    end
    tid{t} = repmat(ismember(t,[1,2,6,7,26,27,28]),[numel(units{t}),1]);
end

rwpa = cat(1,rweightedPos{:});
tida = cat(1,tid{:});

p = 4;
figure();
    for a = 1:5,
        subplot2(5,1,a,1);
        hold('on');
        for i = find(~tida)',
            plot(sq(rwpa(~tida,2,[1,p],a))',sq(rwpa(~tida,1,[1,p],a))');
        end
    end



figure();
for p = 1:5;
    for a = 1:5,
        subplot2(5,5,a,p);
        hold('on');
        plot(sq(rwpa(~tida,2,p,a))',sq(rwpa(~tida,1,p,a))','.');
        plot(mean(sq(rwpa(~tida,2,p,a))'),mean(sq(rwpa(~tida,1,p,a))'),'or');
        xlim([-150,150]);
        ylim([-150,150]);        
        grid('on');
    end
end    
        

    
figure();
%hba < -0.4
% IN   bodyMazeAng.data > -pi/4  & bodyMazeAng.data < pi/4;
% OUT (bodyMazeAng.data > pi*3/4 | bodyMazeAng.data < -pi*3/4);
% CW  -pi*3/4 < bodyMazeAng.data & bodyMazeAng.data < -pi/4;
% CCW  pi*1/4 < bodyMazeAng.data & bodyMazeAng.data < pi*3/4;
%ndgrid(pfl{t}.adata.bins{1:2})
t = 24; 
uint = 0;
for u = 1:numel(units{t}),
    uint = uint+1;
    clf();    
    rmax = max(cell2mat(cf(@(p) prctile(p.data.rateMap(:,p.data.clu==units{t}(u),1),99.5), pfl{t})));
    for p = 1:numel(pfl{t}),
        rmap = plot(pfl{t}{p},units{t}(u));
        al = 1:5;
        for a = 1:numel(al)
            
            sax(a,p) = ...
            subplot2(numel(al),numel(pfl{t}),a,p);
                pcolor(pfl{t}{p}.adata.bins{1},...
                       pfl{t}{p}.adata.bins{2},...
                       rmap(:,:,al(a)));
                caxis([0,rmax]);
                colormap('jet');
                shading('flat');
                axis('xy');
                Lines([],0,'k');
                Lines(0,[],'k');
                xlim(lims{1});
                ylim(lims{2});        
                if p ==1,    
                    ylabel(num2str(pfl{t}{p}.adata.bins{3}(al(a))));
                end
        end
    end
    title(num2str(units{t}(u)));
    waitforbuttonpress();
end

figure,
    hold('on');
    imagesc(pfl{t}{p}.adata.bins{1},...
            pfl{t}{p}.adata.bins{2},...
            rmap(:,:,a)');
    plot3(rweightedPos(1,1,p,a),rweightedPos(1,2,p,a),10,'*m');
    

figure
 data = rand(10);
 data(data > 0.9) = NaN;
 [nr,nc] = size(data);
 subplot(2,1,1); 
 imagesc(data); 
 subplot(2,1,2); 
 pcolor([data nan(nr,1); nan(1,nc+1)]);

 set(gca, 'ydir', 'reverse');

% $$$ compute_egohba_ratemap(Trial,units,xyz,spk,pft,thetaRefChan)
% sessionList, Trials, units, xyz, spk, pft
% $$$ for tind = 1:23; % jg05-20120312.cof.all
% $$$ 
% $$$     %function [pfs] = compute_egohba_pfs(Trial,units,xyz,spk,pft)
% $$$ 
% $$$ 
% $$$     if isempty(units),
% $$$         continue;
% $$$     end;
% $$$ 
% $$$     sampleRate = xyz.sampleRate;
% $$$ 
% $$$ % COMPUTE anglular difference between the head and body
% $$$     headBodyAng = [xyz(:,'spine_upper',[1,2])-xyz(:,'bcom',[1,2]),...
% $$$                    xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2])];
% $$$     headBodyAng = sq(bsxfun(@rdivide,headBodyAng,sqrt(sum(headBodyAng.^2,3))));
% $$$     headBodyAng = cart2pol(headBodyAng(:,:,1),headBodyAng(:,:,2));
% $$$     headBodyAng = circ_dist(headBodyAng(:,2),headBodyAng(:,1));
% $$$     headBodyAng = MTADfet.encapsulate(Trial,...
% $$$                                       -(headBodyAng+hbaCorrection(tind)),...
% $$$                                       sampleRate,...
% $$$                                       'hba','hba','h');
% $$$     
% $$$ % TRANSFORM Local Field Potential -> theta phase
% $$$     phz = load(Trial,'lfp',thetaRefChan).phase([6,12]);
% $$$     phz.data = unwrap(phz.data);
% $$$     phz.resample(xyz);    
% $$$     phz.data = mod(phz.data+pi,2*pi)-pi;
% $$$ 
% $$$     hvec = xyz(:,'nose',[1,2])-xyz(:,'hcom',[1,2]);
% $$$     hvec = sq(bsxfun(@rdivide,hvec,sqrt(sum(hvec.^2,3))));
% $$$     hvec = cat(3,hvec,sq(hvec)*[0,-1;1,0]);
% $$$     hvec = multiprod(hvec,...
% $$$                      [cos(rot(t)),-sin(rot(t));sin(rot(t)),cos(rot(t))],...
% $$$                      [2,3],...
% $$$                      [1,2]);
% $$$     
% $$$     
% $$$     thetaState = resample(cast([Trial.stc{'theta-groom-sit-rear'}],'TimeSeries'),xyz);
% $$$     
% $$$     pfTemp = Trial;
% $$$ 
% $$$     pargs = get_default_args('MjgER2016','MTAApfs','struct');        
% $$$     pargs.units        = units;
% $$$     pargs.tag          = 'egofield';
% $$$     pargs.binDims      = [20, 20, 0.1];
% $$$     pargs.SmoothingWeights = [1.8, 1.8, 1.8];
% $$$     pargs.halfsample   = false;
% $$$     pargs.numIter      = 1;   
% $$$     pargs.boundaryLimits = [-400,400;-400,400;-pi,pi];
% $$$     pargs.states       = '';
% $$$     pargs.overwrite    = true;
% $$$     pargs.autoSaveFlag = false;    
% $$$     electrode = 0;
% $$$     
% $$$     for p = 1:numel(binPhzc)
% $$$         pargs.tag          = ['egofield_theta_phase_',num2str(p)];
% $$$         for u = 1:numel(units),
% $$$             if u==1 | electrode~=spk.map(spk.map(:,1)==units(u),2), % update phase state            
% $$$                 pargs.spk = copy(spk);
% $$$                 electrode = spk.map(spk.map(:,1)==units(u),2);
% $$$                 pargs.states = copy(thetaState);
% $$$                 pargs.states.label = ['thetaPhz_',num2str(p)];
% $$$                 pargs.states.data(   (phz(:,electrode) < binPhzs(p)   )    ...
% $$$                                    | (phz(:,electrode) >= binPhzs(p+1)) ) = 0;
% $$$                 cast(pargs.states,'TimePeriods');
% $$$                 resInd = WithinRanges(pargs.spk.res,pargs.states.data);
% $$$                 pargs.spk.res = pargs.spk.res(resInd);
% $$$                 pargs.spk.clu = pargs.spk.clu(resInd);
% $$$             end 
% $$$             
% $$$             [mxr,mxp] = pft.maxRate(units(u));
% $$$             pfsCenterHR = MTADfet.encapsulate(Trial,                                           ...
% $$$                                               [multiprod(bsxfun(@minus,                        ...
% $$$                                                                 mxp,                           ...
% $$$                                                                 sq(xyz(:,'hcom',[1,2]))),      ...
% $$$                                                          hvec,2,[2,3]),headBodyAng],           ...
% $$$                                               sampleRate,                                      ...
% $$$                                               'egocentric_placefield',                         ...
% $$$                                               'egopfs',                                        ...
% $$$                                               'p'                                              ...
% $$$                                               );
% $$$             pargs.xyzp = pfsCenterHR;
% $$$             pargs.units  = units(u);
% $$$             pfsArgs = struct2varargin(pargs);
% $$$             pfTemp = MTAApfs(pfTemp,pfsArgs{:});
% $$$             if u==1,
% $$$                 pfTemp.purge_savefile();
% $$$                 pfTemp.save();        
% $$$             end    
% $$$         end
% $$$         pfTemp.save();
% $$$         pfTemp = Trial;
% $$$     end
% $$$ end
% $$$ 
% $$$ %func end

tind = 18;
Trial = Trials{tind};
unitSubset = units{tind};
tprpf = {};
for p = 1:numel(binPhzc),
    tprpf{p} = MTAApfs(Trial,[],[],[],['egofield_theta_phase_',num2str(p)]);
end
pft = pfs_2d_theta(Trial,unitSubset);

figure,
hax = tight_subplot(1,numel(binPhzc)+1,0.01,0.2,0.1);
for u = unitSubset,
    for p = 1:numel(binPhzc),            
        axes(hax(p));
        cla();
        plot(tprpf{p},u,1,'',[0,15],false);
        title({['unit: ' num2str(u)],['phzBin: ' num2str((binPhzc(p))/pi*180)]});        
    end
    axes(hax(end));
    h = gca();
    cla();
    plot(pft,u,'mean','colorbar',[],true);    
    h.Position(3) = hax(end-1).Position(3);
    waitforbuttonpress();
end


tprpf = {};
for tind = 17:23;
    Trial = Trials{tind};
    unitSubset = units{tind};
    for p = 1:numel(binPhzc),
        tprpf{tind-16,p} = MTAApfs(Trial,[],[],[],['egofield_theta_phase_',num2str(p)]);
    end
end


mpfsr = repmat({[]},[1,numel(binPhzc)]);
for tind = 1:size(tprpf,1),
    for p = 1:numel(binPhzc),
    mpfsr{p} = cat(2,mpfsr{p},tprpf{tind,p}.data.rateMap);
    end
end


cf(@(t) t.load('nq'), Trials(17:23));
edist = cf(@(t,u) t.nq.eDist(u), Trials(17:23),units(17:23));
edist = cat(1,edist{:});


mpfsrm = cf(@(m,p,e) reshape(mean(m(:,max(m)'>1&e>20),2,'omitnan'),p.adata.binSizes'),mpfsr,tprpf(1,:),repmat({edist},[1,numel(mpfsr)]));

figure,

clf();
hax = tight_subplot(1,numel(binPhzc),0.01,0.2,0.1);
for p = 1:numel(binPhzc),            
    axes(hax(p));
    cla();
    imagesc(tprpf{p}.adata.bins{:},mpfsrm{p}');
    title({['phzBin: ' num2str((binPhzc(p))/pi*180)]});        
    caxis([0,7.5]);
    grid('on');
    clear_axes_labels(gca);    
    lax = Lines(0,[],'r');
    if p == 1,
        hax(p).XAxisLocation = 'bottom';        
        xlabel('longitudinal (back/front)');
        ylabel('lateral (left/right)');
    end
    axis('xy');    
end
ForAllSubplots('h = gca();h.Units=''centimeters'';h.Position(4) = h.Position(3);');
ForAllSubplots('h = gca();h.GridColor=[1,1,1];');

fax = axes();
fax.Position = [0,0,1,1];
fax.Visible = 'off';
xlim([0,1]);
ylim([0,1]);
tax = text(fax,0.5,0.95,'Mean firing rate of place field in egocentric coordinates partitioned by theta phase. Subject:jg05','HorizontalAlignment','center');

hax(end).Units = 'normalized';
line(hax(end).Position(1)+hax(end).Position(3).*[0.75,1],...
     [0.1,0.1],'LineWidth',1,'Color',[0,0,0]);
text(hax(end).Position(1)+hax(end).Position(3).*0.75,...
     0.05,'10 cm');
cax = colorbar(hax(end));
cax.Position(1) = cax.Position(1)+0.03;

print(gcf,'-dpng',fullfile('/storage/share/Projects/BehaviorPlaceCode/phase_precession',...
                             'ego_pp_meanFiringRate.png'));
print(gcf,'-depsc2',fullfile('/storage/share/Projects/BehaviorPlaceCode/phase_precession',...
                             'ego_pp_meanFiringRate.eps'));




%        PROCESS : ESI <- information <- rate maps <- MTAApfs objects

% MTAApfs objects: pfsh pfs
pfs1 = cf(@(p) p{4},pfs);

bins = pfs1{1}.adata.bins;
binSizes = pfs1{1}.adata.binSizes;
width = binSizes(1);
height =binSizes(2);
radius = round(binSizes(1)/2)-find(bins{1}<=-250,1,'last');
centerW = width/2;
centerH = height/2;
[W,H] = meshgrid(1:width,1:height);           
mazeMask = logical(reshape(double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius),[],1));

[erm,clu] = decapsulate_and_concatenate_mtaapfs(pfs1,units);
erm = reshape(erm,[prod(binSizes(1:2)),binSizes(3),size(clu,1)]);
erm = erm(mazeMask,:,:);


esi = sq(sum(1/size(erm,1).*bsxfun(@rdivide,erm,mean(erm)) ...
          .*log2(bsxfun(@rdivide,erm,mean(erm))),'omitnan'))';


pfsh1 = cf(@(p) p{4},pfsh);
ermSH = decapsulate_and_concatenate_mtaapfs(pfsh1,units);
ermSH = reshape(ermSH,[prod(binSizes(1:2)),binSizes(3),size(clu,1),size(ermSH,3)]);
ermSH = ermSH(mazeMask,:,:,:);

esiSH = permute(sq(sum(1/size(ermSH,1).*bsxfun(@rdivide,ermSH,mean(ermSH)) ...
          .*log2(bsxfun(@rdivide,ermSH,mean(ermSH))),'omitnan')),[2,1,3]);


esiZ = (log2(esi) - mean(log2(esiSH),3))./std(log2(esiSH),[],3);
u = 149;
figure();
hist(log2(sq(esiSH(u,2,:))),20);
Lines(log2(esi(u,2)),[],'r');


figure
u = 25;
for a = 1:5
    subplot2(2,5,1,a);
    plot(pfs{20}{4},u,{':',':',a,1},'text',[],false,[]);    
    Lines([],0,'k');
    xlim([-250,250]);
    ylim([-250,250]); 
    subplot2(2,5,2,a);
    plot(pfsh{20}{4},u,{':',':',a,3},'text',[],false,[]);
    Lines([],0,'k');    
    xlim([-250,250]);
    ylim([-250,250]);    
end


 
% GENERATE computational mask 
lims = {[-250,250],[-250,250]};
maskBinInd = cf(@(p,l) l(1) < p & p < l(2),  pfs{t}{1}.adata.bins([1:2]),lims);
mask = maskBinInd{1};
for dim = 2:numel(maskBinInd),
     mask = bsxfun(@and,mask,permute(maskBinInd{dim},[2:dim,1]));
end
maskDimVal = cf(@(p,i) p(i), pfs{t}{1}.adata.bins(1:2),maskBinInd);
maskBinPos = cell(size(maskDimVal));
[maskBinPos{:}] = ndgrid(maskDimVal{:});
maskBinPos = reshape(cat(ndims(maskBinPos)+1,maskBinPos{:}),[],ndims(maskBinPos));

rateMapNormal = pfs;
rateMapShuffled = pfsh;
rateMapNormal = pfl;
rateMapShuffled = pflh;


% COMPUTE the rate weighted field postion
rweightedPos = {};
trialAnatomicalGroup = {};

goodTrialInds = find(~cellfun(@isempty,units));

nAng = 5;
for trialId = goodTrialInds,
    rweightedPos{trialId} =                                                               ...
        nan([numel(units{trialId}),                                                       ...
             2,                                                                           ...
             numel(rateMapNormal{trialId}),                                                     ...
             rateMapNormal{trialId}{1}.adata.binSizes(end)]);
    
    for unitId = 1:numel(units{trialId}),
        for thetaPhase = 1:numel(rateMapNormal{trialId}),
            rmap = plot(rateMapNormal{trialId}{thetaPhase},units{trialId}(unitId));
            rmap(rmap<1)=0;
            % Ratemap weighted position of egocetric placefield
            %     maskBinPos - coordinates of each bin
            %     ratemap normalized to 1
            for hbAngle = 1:nAng
                rweightedPos{trialId}(unitId,:,thetaPhase,hbAngle) =                      ...
                    sum(bsxfun(@times,                                                    ...
                               maskBinPos,                                                ...                           
                               reshape(rmap(maskBinInd{:},hbAngle),[],1)                  ...
                               ./sum(reshape(rmap(maskBinInd{:},hbAngle),[],1),'omitnan'))...
                        ,'omitnan');
            end%for hbAngle
        end%for thetaPhase
    end%for unitId
    trialAnatomicalGroup{trialId} = repmat(ismember(trialId,[1,2,6,7,26,27,28]),[numel(units{trialId}),1]);
end%for trialId

% rwpa
rateWeightedPosition = cat(1,rweightedPos{:});
trialAnatomicalGroup = cat(1,trialAnatomicalGroup{:});

rweightedPosSH = {};
nAng  = 5;
nIter = size(rateMapShuffled{1}{1}.data.rateMap,3);
for trialId = goodTrialInds,
    rweightedPosSH{trialId} =                                                             ...
        nan([numel(units{trialId}),                                                       ...
             2,                                                                           ...
             numel(rateMapShuffled{trialId}),                                                   ...
             rateMapShuffled{trialId}{1}.adata.binSizes(end),                                   ...
             nIter]);
    for unitId = 1:numel(units{trialId}),
        for thetaPhase = 1:numel(rateMapShuffled{trialId}),
            clusterIndex = units{trialId}(unitId)==pfsh{trialId}{thetaPhase}.data.clu;
            rmap = reshape(pfsh{trialId}{thetaPhase}                                      ... MTAApfs object
                             .data                                                        ... data
                               .rateMap(:,                                                ... rate values
                                        clusterIndex,                                     ... cluster Index
                                        :),                                               ... samples
                           [],nAng,nIter);
            rmap(rmap<1)=0;            
            % Ratemap weighted position of egocetric placefield
            %     maskBinPos - coordinates of each bin
            %     ratemap normalized to 1
            rweightedPosSH{trialId}(unitId,:,thetaPhase,:,:) =                            ...
                permute(sum(bsxfun(@times,                                                ...
                                repmat(maskBinPos,[1,1,nAng,nIter]),                      ...
                                permute(rmap(mask(:),:,:),[1,4,2,3])                      ...
                                   ./ sum(permute(rmap(mask(:),:,:),[1,4,2,3]),           ...
                                          'omitnan')),                                    ...
                            'omitnan'),                                                   ...
                        [1,2,5,3,4]);
        end%for thetaPhase
    end%for unitId
end%for trialId



unitSessionMap = cf(@(u,i) u*0+i, units,num2cell(1:numel(units)));
unitSessionMap = cat(2,cat(2,unitSessionMap{:})',cat(2,units{:})');

%rwpaSH = cat(1,rweightedPosSH{:});
rateWeightedPositionShuffled = cat(1,rweightedPosSH{:});

zscr = (bsxfun(@minus,                                                                    ...
               rateWeightedPosition,                                                      ...
               mean(rateWeightedPositionShuffled,5)))                                     ...
       ./std(rateWeightedPositionShuffled,[],5);


zscrFWE = bsxfun(@minus,                                                                    ...
               rateWeightedPositionShuffled,                                              ...
               repmat(mean(rateWeightedPositionShuffled,5),[1,1,1,1,size(rateWeightedPositionShuffled,5)]))                                     ...
       ./repmat(std(rateWeightedPositionShuffled,[],5),[1,1,1,1,size(rateWeightedPositionShuffled,5)]);


p  = 1;
x  = 2;
al = 2;
ar = 4;
tida   = logical(trialAnatomicalGroup);
rwpa   = rateWeightedPosition;
rwpaSH = rateWeightedPositionShuffled;

figure();
for a = 1:nAng,
    subplot2(5,1,a,1);
    hist(zscr(~trialAnatomicalGroup,2,4,a),15);
    xlim([-6,6]);
end

figure();
for a = 1:nAng,
    subplot2(5,1,a,1);
    hist(rwpa(~trialAnatomicalGroup,2,4,a),15);
    xlim([-120,120]);
end




sigThresh = 1.94;
sum(zscr(~tida,x,p,ar)>sigThresh & zscr(~tida,x,p,al)<-sigThresh)

p2z = @(p) icdf('normal',p,0,1);



plot(zscrFWE(~tida,x,p,al),zscrFWE(~tida,x,p,ar),'.m');    

plot(zscrFWE(~tida,x,p,al),zscrFWE(~tida,x,p,ar),'.m');    

rotation_matrix = @(t) [cos(t),-sin(t);sin(t),cos(t)];

vxfrmZ = plot(zscr(~tida,x,p,al),zscr(~tida,x,p,ar),'.');
xfrmZFWE = (rotation_matrix(pi/4)*sq(zscrFWE(~tida,x,p,[al,ar],1))')';


xfrmT = p2z(0.95)*std(xfrmZFWE(:,1))+mean(xfrmZFWE(:,1));
figure();
plot(xfrmZFWE(:,1),xfrmZFWE(:,2),'.');

pc = [];
r = [];
for i = 1:100,
    [R_CA1, Pval_CA1] = corrcoef(zscrFWE(~tida,x,p,al,i),zscrFWE(~tida,x,p,ar,i));
    [P_CA1,S_CA1]     = polyfit (zscrFWE(~tida,x,p,al,i),zscrFWE(~tida,x,p,ar,i),1);
    pc(i,:) = P_CA1;
    r(i) = R_CA1(2);
end;%for i
figure();
hist(r,20);
Lines(R_CA1(2),[],'r');

figure();
hist(pc(:,1),20);
Lines(P_CA1(1),[],'r');


figure();
hold('on');
plot(pc(:,1),pc(:,2),'.');
[P_CA1,S_CA1] = polyfit(zscr(~tida,x,p,al),zscr(~tida,x,p,ar),1);    
[R_CA1, Pval_CA1] = corrcoef(zscr(~tida,x,p,al),zscr(~tida,x,p,ar));
plot(P_CA1(:,1),P_CA1(:,2),'*r');

x = 2;
ny = 1;
figure();
%for p = 1:5,
for p = 4,
    y = 1;
    subplot2(ny,2,y,1);
    hold('on');
    grid('on');    
    % CA1
    sigThresh = p2z(1-(1-0.95)^(1/sum(~tida)));
    plot(zscr(~tida,x,p,al),zscr(~tida,x,p,ar),'.b');
    plot(zscrFWE(~tida,x,p,al,1),zscrFWE(~tida,x,p,ar,1),'.m');    
    plot(mean(zscr(~tida,x,p,al)),mean(zscr(~tida,x,p,ar)),'ob');
    [P_CA1,S_CA1] = polyfit(zscr(~tida,x,p,al),zscr(~tida,x,p,ar),1);
    line([-6,6],polyval(P_CA1,[-6,6]),'Color','b');
    [R_CA1, Pval_CA1] = corrcoef(zscr(~tida,x,p,al),zscr(~tida,x,p,ar));
    % Subplot Opts    
    xlim([-6,6]);
    ylim([-6,6]);
    Lines([],-sigThresh,'k','--','LineWidth',1);
    Lines(sigThresh,[],'k','--','LineWidth',1);
    line([-6,6],[6,-6],'Color','c');
    xlabel('leftward HBA (z-score)')
    ylabel('rightward HBA (z-score)')
    title('CA1 ');
    
    subplot2(ny,2,y,2);
    hold('on');
    grid('on');        
    % CA3
    sigThresh = p2z(1-(1-0.95)^(1/sum(tida)));
    plot(zscr(tida,x,p,al),zscr(tida,x,p,ar),'.g');
    plot(zscrFWE(tida,x,p,al,1),zscrFWE(tida,x,p,ar,1),'.m');        
    plot(mean(zscr(tida,x,p,al)),mean(zscr(tida,x,p,ar)),'og');
    [P_CA3, S_CA3   ] = polyfit(zscr(tida,x,p,al),zscr(tida,x,p,ar),1);
    line([-6,6],polyval(P,[-6,6]),'Color','g');
    [R_CA3, Pval_CA3] = corrcoef(zscr(tida,x,p,al),zscr(tida,x,p,ar));
    % Subplot Opts        
    xlim([-6,6]);
    ylim([-6,6]);
    Lines([],-sigThresh,'k','--','LineWidth',1);
    Lines(sigThresh,[],'k','--','LineWidth',1);
    line([-6,6],[6,-6],'Color','c');
    xlabel('leftward HBA (z-score)')
    ylabel('rightward HBA (z-score)')
    title('CA3 ');
end


figure();
p = 4;
for a = 1:nAng, for x = 1:2,
subplot2(nAng,4,a,x);
    hold('on');
    mrwpa(x,a,p) = mean(rwpa(~tida,x,p,a));
    plot(rwpa(~tida,x,p,a),zscr(~tida,x,p,a),'.b');
    plot(rwpa(tida,x,p,a), zscr(tida,x,p,a), '.r');
    xlabel(['mean pos (cm)']);
    ylabel(['z-score (A.U.)']);    
    grid('on');
    if (x==2),
        legend({'CA1','CA3'},'location','southeast');
    end
    %title();
    xlim([-150,150]);
    ylim([-6,6]);    
end;end;
subplot2(nAng,4,[1:2],[3:4]);


sq(sum(zscr(tida,x,:,ar)>2&-2>zscr(tida,x,:,al)))./sum(tida)*100

sq(sum(zscr(~tida,x,:,ar)>2|-2>zscr(~tida,x,:,al)))./sum(~tida)*100

zu = find( (-2> zscr(:,x,p,al)) & (zscr(:,x,p,ar) > 2) );
 
figure();
%for z = 60:170,%zu'    
for z = zu'

t = clu(z,1);
u = find(units{t}==clu(z,2));
pA = 4;
pT = 3;
pD = 2;
rmapA = plot(pfl{t}{pA},units{t}(u));
rmapD = plot(pfl{t}{pD},units{t}(u));
rmapT = plot(pfl{t}{pT},units{t}(u));
for a = 1:nAng
    subplot2(5,4,a,1);
        hold('on');
        hpc = pcolor(pfs{t}{pA}.adata.bins{1:2},rmapD(:,:,a));
        hpc.EdgeColor = 'none';
        plot(sq(rwpaSH(z,2,pD,a,:)),sq(rwpaSH(z,1,pD,a,:)),'.c');
        plot(sq(  rwpa(z,2,pD,a,:)),sq(  rwpa(z,1,pD,a,:)),'*r');
        xlim([-400,400]);
        ylim([-400,400]);
        caxis([0,12]);
    subplot2(5,4,a,2);
        hold('on');
        hpc = pcolor(pfs{t}{pT}.adata.bins{1:2},rmapT(:,:,a));
        hpc.EdgeColor = 'none';
        plot(sq(rwpaSH(z,2,pT,a,:)),sq(rwpaSH(z,1,pT,a,:)),'.c');
        plot(sq(  rwpa(z,2,pT,a,:)),sq(  rwpa(z,1,pT,a,:)),'*r');
        xlim([-400,400]);
        ylim([-400,400]);
        caxis([0,12]);
    subplot2(5,4,a,3);
        hold('on');
        hpc = pcolor(pfs{t}{pA}.adata.bins{1:2},rmapA(:,:,a));
        hpc.EdgeColor = 'none';
        plot(sq(rwpaSH(z,2,pA,a,:)),sq(rwpaSH(z,1,pA,a,:)),'.c');
        plot(sq(  rwpa(z,2,pA,a,:)),sq(  rwpa(z,1,pA,a,:)),'*r');
        xlim([-400,400]);
        ylim([-400,400]);
        caxis([0,12]);
end
subplot2(5,4,1,4);
plot(pft{t},units{t}(u),1,'text',[0,6]);
subplot2(5,4,2,4);
plot(pfe{t},units{t}(u),1,'text',[0,6],false);
title(num2str([z,clu(z,:)]));
waitforbuttonpress();
end

17/111

1-(1-0.95)^(1/111)



figure();
for u = 80
for a = 1:nAng,
for p = 1:nAng,    
    subplot2(5,5,p,a);
    hold('on');
    plot(sq(rwpaSH(u,2,4,a,:)),sq(rwpaSH(u,1,4,a,:)),'.');
    plot(sq(  rwpa(u,2,4,a,:)),sq(  rwpa(u,1,4,a,:)),'*r');
    xlim([-200,200]);
    ylim([-200,200]);
    grid('on');
end
title(num2str([u,clu(u,:)]));
waitforbuttonpress();
clf();
end


%%%<< plot ego ratemap
figure()
t = 20;
for u = 1:numel(units{t});
unit = units{t}(u);
pcolor(pfe{t}.adata.bins{2},...
       pfe{t}.adata.bins{1},...
       plot(pfe{t},unit,[],[],[],false));
        
colormap(gca(),'jet');
shading (gca(),'flat');
axis    (gca(),'xy');
xlim    (gca(),lims{1});
ylim    (gca(),lims{2});        

Lines([],0,'k');
Lines(0,[],'k');
set(gca(),'XTick',[]);
set(gca(),'YTick',[]);

    
title(['egoratemap unit:',num2str(unit)]);
waitforbuttonpress();
clf();
end
%%%>>>


%%%<< plot ego theta ratemaps
figure()
t = 4;
for u = 1:numel(units{t});
    unit = units{t}(u);
    for p = 1:5,
        subplot(1,5,p)
        pcolor(pfet{t}{p}.adata.bins{2},...
               pfet{t}{p}.adata.bins{1},...
               plot(pfet{t}{p},unit,[],[],[],false));
        colormap(gca(),'jet');
        shading (gca(),'flat');
        axis    (gca(),'xy');
        xlim    (gca(),lims{1});
        ylim    (gca(),lims{2});        
        Lines([],0,'k');
        Lines(0,[],'k');
        set(gca(),'XTick',[]);
        set(gca(),'YTick',[]);
    end    
    title(['egoratemap unit:',num2str(unit)]);
    waitforbuttonpress();
    clf();
end
%%%>>>
