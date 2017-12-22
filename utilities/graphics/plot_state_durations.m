function plot_state_durations(Trial,varargin)

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('stcMode',  'msnn_ppsvd_raux'                                                   ...
);
[stcMode] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

stc = Trial.load('stc',stcMode);


%% First Level
statesPrime = {'loc','rear','pause','groom', 'sit' };
sc          = [0,0,1; 1,0,0;  0,1,1;  1,0,1;  1,1,0]*.5;  %yellow
% COMPUTE state occupancy
spdur = cell2mat(cf(@(s,sts) sum(diff([s{sts,1}.data],1,2)),repmat({stc},[1,numel(statesPrime)]),statesPrime));
% COMPUTE states' positions as blocks ordered left to right
sppos = [0,cumsum(spdur(1:end-1))];
% PLOT Primary Level: all states mutually exclusive 
for s = 1:numel(statesPrime);patch(sppos(s)+[0,spdur(s),spdur(s),0],[2,2,1,1],sc(s,:));end    


%% Second Level
statesSecond = {'loc&theta','rear&theta','pause&theta','groom&theta','sit&theta'};
% DARK for non theta blocks
sc           = [      0,0,1;       1,0,0;        0,1,1;        1,0,1;      1,1,0]*0.4+0.25;  
% LIGHT for theta blocks
sct          = [      0,0,1;       1,0,0;        0,1,1;        1,0,1;      1,1,0];
% COMPUTE State durations
ssdur = cell2mat(cf(@(s,sts) sum(diff([s{sts,1}.data],1,2)),repmat({stc},[1,numel(statesSecond)]),statesSecond));

gdur = reshape([ssdur;spdur-ssdur],[],1);
gpos = [0;cumsum(gdur(1:end-1))];
scg  = sq(reshape(permute(cat(3,sct,sc),[3,1,2]),[],1,3));
% PLOT Secondary Level: all states mutually exclusive 
for s = 1:numel(gdur); patch(gpos(s)+[0,gdur(s),gdur(s),0],[1,1,0,0],scg(s,:));end    

%% Third Level
statesThird = {'lloc',      'hloc','rear','lpause',    'hpause', 'groom', 'sit'};
sc          = [ 0,0,1; .25,.25,.75; 1,1,1;   0,1,1; .25,.75,.75;  1,1,1;  1,1,1];        
% COMPUTE state occupancy
spdur = cell2mat(cf(@(s,sts) sum(diff([s{sts,1}.data],1,2)),repmat({stc},[1,numel(statesThird)]),statesThird));
% COMPUTE states' positions as blocks ordered left to right
sppos = [0,cumsum(spdur(1:end-1))];
% PLOT Primary Level: all states mutually exclusive 
for s = 1:numel(spdur); patch(sppos(s)+[0,spdur(s),spdur(s),0],[0,0,-1,-1],sc(s,:)); end    


%% Fourth Level
statesFourth = {'lloc&theta','lloc-theta','hloc&theta','hloc-theta',...
               'rear','lpause&theta','lpause-theta','hpause&theta','hpause-theta',...
               'groom','sit'};
sc  = [0,0,1;         ... l blue   lloc&theta
       0.5,0.5,1;     ... l violet lloc-theta
       0.25,0.25,0.75;... l cyan 
       0.5,0.5,1;     ... l violet lloc-theta
       1,1,1;         ... white    
       0,1,1;         ... cyan
       0,0.5,0.5;     ... dark cyan              
       0.25,0.75,0.75;... light cyan
       0,0.5,0.5;     ... dark cyan       
       1,1,1;         ... white
       1,1,1];          % white
% COMPUTE state occupancy
spdur = cell2mat(cf(@(s,sts) sum(diff([s{sts,1}.data],1,2)),repmat({stc},[1,numel(statesFourth)]),statesFourth));
% COMPUTE states' positions as blocks ordered left to right
sppos = [0,cumsum(spdur(1:end-1))];

% PLOT Primary Level: all states mutually exclusive 
for s = 1:numel(spdur);
    patch(sppos(s)+[0,spdur(s),spdur(s),0],... X
          [-1,-1,-2,-2],      ... Y
          sc(s,:));         % C
end    


set(gca,'YTick',[-1.5,-0.5,0.5,1.5]);
set(gca,'YTickLabels',{'HL theta','High/Low','theta','State'});
axis('tight');
set(gca,'XTick',[]);