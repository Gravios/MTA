% MjgER2016 Figure3
%
% Subplots:
%    A. autoccgs
%    B. Placefield examples for each state
%    C. 
%    D. Feature Matrix
%    E. Behavioral labels - hand labels and neural network labels



sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';

FigDir = '/storage/gravio/figures/placefields';
mkdir(FigDir);

% SET behvaioral states
states = {'loc&theta','lloc&theta','hloc&theta','rear&theta',    ...
          'pause&theta','lpause&theta','hpause&theta',           ...
          'theta-groom-sit'};
statesCcg = {'loc','lloc','hloc','rear',    ...
          'pause','lpause','hpause',           ...
          'theta-groom-sit'};

 
% LOAD Trials
% COMPUTE placefield statistics
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);
pfstats = cf(@(t)  compute_pfstats_bs(t),       Trials);


% $$$ cf(@(t)  MjgER2016_drzfields(t,true), Trials);




patchComparison = [];

% FOR each Trial -------------------------------------------------------------
for t = 1:20,

% LOAD Trial
    Trial = Trials{t};
    units = pfstats{t}.cluMap;    
    disp(['Processing trial: ' Trial.filebase]);

% LOAD placefields and subsampled estimate
    for sts = 1:numel(states),        
        defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
        defargs.units = units;
        %defargs.overwrite = true;
        defargs.states = states{sts};
        defargs = struct2varargin(defargs);        
        pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
    end

% COMPARE patches
    for u = 1:numel(units),
        unit = units(u);
        for s = 1:8,
            for o = 1:8,
                rmap = plot(pfkbs{o},unit,'mean');
                patchInd = sq(pfstats{t}.pfmstats(s,unit==units).patchRateInd(1,1,1,:,:));
                oPatchRates = rmap(sub2ind(size(rmap),patchInd(1,nniz(patchInd(1,:)'))',...
                                           patchInd(2,nniz(patchInd(2,:)'))'));
                if s == 1 && o ==1, pcuInd = size(patchComparison,1)+1; end
                
                if sum(nniz(oPatchRates))>0.9*numel(oPatchRates), 
                    patchComparison(pcuInd,s,o) = mean(oPatchRates,'omitnan');
                    patchComparison(pcuInd,s,o) = mean(oPatchRates,'omitnan');                    
                else,
                    patchComparison(pcuInd,s,o) = nan;
                end
                    
            end%for o
        end%for s 
    end%for u
end%for t

        
figure();
subplot(141); plot(patchComparison(:,4,1),patchComparison(:,4,4),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
xlabel('mean loc rate in rear');ylabel('mean rear rate in rear');
subplot(142); plot(patchComparison(:,1,1),patchComparison(:,1,4),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
subplot(143); plot(patchComparison(:,2,2),patchComparison(:,2,3),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
subplot(144); plot(patchComparison(:,3,2),patchComparison(:,3,3),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');

figure();
subplot(121); plot(patchComparison(:,8,1),patchComparison(:,8,4),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
subplot(122); plot(patchComparison(:,8,2),patchComparison(:,8,3),'.'); daspect([1,1,1]);xlim([0,35]);ylim([0,35]);grid('on');
title;xlabel('low loc');ylabel('high loc');


figure();
subplot(
hist((patchComparison(:,8,2)-patchComparison(:,8,3))./(patchComparison(:,8,2)+patchComparison(:,8,3)),50);
        
% $$$ 
% $$$     
% $$$ 
% $$$     
% $$$     
% $$$ % LOAD DRZ fields
% $$$     dfs = cell([1,3]);
% $$$     [dfs{:}] = MjgER2016_drzfields(Trial,true);
% $$$     dfst = {'pitch','height','rhm'};
% $$$ 
% $$$ 
% $$$ 
% $$$     %[accg,tbins] = autoccg(Trial);
% $$$ 
% $$$     % PLACEFIELD figure now located in placefield_summary.m
% $$$ end
% $$$ 
% $$$ pr = reshape(cell2mat(cf(@(p) [p.pfmstats.peakFR],pfstats)),8,[]);
% $$$ 
% $$$ 
% $$$ figure();
% $$$ 
% $$$ clf();
% $$$ % loc x rear
% $$$ subplot(221);  hold('on');  plot(pr(1,:),pr(4,:),'.');  daspect([1,1,1]);
% $$$                line([0,40],[0,40]);  line([0,40],[0,20]); line([0,20],[0,40]);
% $$$ subplot(222);  hold('on');  plot(pr(2,:),pr(4,:),'.');   daspect([1,1,1]);
% $$$                line([0,40],[0,40]);  line([0,40],[0,20]); line([0,20],[0,40]);
% $$$                
% $$$ subplot(224);  hold('on');  plot(pr(3,:),pr(4,:),'.');   daspect([1,1,1]);
% $$$                line([0,40],[0,40]);  line([0,40],[0,20]); line([0,20],[0,40]);               