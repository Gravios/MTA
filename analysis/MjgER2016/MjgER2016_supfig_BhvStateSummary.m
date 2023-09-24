% req20171107 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: plot_state_durations.m
%  Description: State occupancy breakdown
%  Bugs: NA


% testing phase
% Trial, stcMode
configure_default_args();
MjgER2016_load_data();

statesPrime = {'walk','rear','turn','pause','groom','sit'};

spdur = {};
for tind = 1:numel(Trials)
    Trial = Trials{tind};
    stc = Trial.load('stc','msnn_ppsvd_raux');

% Theta state occupation:
%     Two bars:
%         Top Bar: Patch with length equal to the total duration spent in theta state
%         Bottom Bar: Patches for each total duration of states within theta state


    % COMPUTE state occupancy
    spdur{tind} = cell2mat(cf(@(s,sts) sum(diff([s{sts,1}.data],1,2)),repmat({stc},[1,numel(statesPrime)]),statesPrime));
    % COMPUTE states' positions as blocks ordered left to right
    tdur(tind) = sum(spdur{tind});
end


[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','portrait',[],1.5,1.5,0.05,0.05);
hfig.PaperPositionMode = 'auto';
[~,sind] = sort(tdur);
sax = gobjects([0,1]);
for tind = 1:numel(Trials)
    sax(tind) = axes('Units','centimeters',...
                     'Position',[2.5,tind/2+8.5,3.5*(tdur(sind(tind))/60/60),0.5]);
spos = [0,cumsum(spdur{sind(tind)}(1:end-1))];
sc  = [0,0,1;... blue
       1,0,0;... red
       0,1,0;... green
       0,1,1;... cyan
       1,0,1;... magenta
       1,1,0];  %yellow
% PLOT Primary Level: all states mutually exclusive 
    for s = 1:numel(statesPrime);
        patch(spos(s)+[0,spdur{sind(tind)}(s),spdur{sind(tind)}(s),0],... X
              [2,2,1,1],      ... Y
              sc(s,:));         % C
    end    
    axis(sax(end),'tight');
    sax(end).XTick = []
    sax(end).YTick = []    
end
axes(fax);
line(sax(end).Position(1).*[1,1]+[0,sax(end).Position(3)],sax(end).Position(2)+1.*[1,1]);
line(sax(end).Position(1).*[1,1]+[0,sax(end).Position(3)./(tdur(sind(tind))/60/60).*4],sax(end).Position(2)+1.*[1,1],'Color','r');
line(sax(end).Position(1).*[1,1]+[0,sax(end).Position(3)./(tdur(sind(tind))/60/60).*2],sax(end).Position(2)+1.*[1,1],'Color','k');


% NEXT part
% State transition probabalities
% feature trace at state transitions
 

statesSecond = {'walk&theta','rear&theta','turn&theta','pause&theta','groom&theta','sit&theta'};

% DARK for non theta blocks
sc  = [0,0,1;... light blue
       1,0,0;... red
       0,1,0;... green
       0,1,1;... cyan
       1,0,1;... magenta
       1,1,0]*0.4+0.25;  %yellow

% LIGHT for theta blocks
sct = [0,0,1;... blue
       1,0,0;... red
       0,1,0;... green
       0,1,1;... cyan
       1,0,1;... magenta
       1,1,0];  %yellow

ssdur = cell2mat(cf(@(s,sts) sum(diff([s{sts,1}.data],1,2)),repmat({stc},[1,numel(statesSecond)]),statesSecond));

gdur = reshape([ssdur;spdur-ssdur],[],1);
gpos = [0;cumsum(gdur(1:end-1))];
scg  = sq(reshape(permute(cat(3,sct,sc),[3,1,2]),[],1,3));

% PLOT Secondary Level: all states mutually exclusive 
for s = 1:numel(gdur);
    patch(gpos(s)+[0,gdur(s),gdur(s),0],... X
          [1,1,0,0],      ... Y
          scg(s,:));         % C
end    


% $$$ 
% $$$ statesThird = {'lloc','hloc','rear','lpause','hpause','groom','sit'};
% $$$ sc  = [0,0,1;... light blue
% $$$        0.25,0.25,0.75; ... blue
% $$$        1,1,1;          ... white
% $$$        0,1,1;      ... cyan
% $$$        0.25,0.75,0.75; ... cyan
% $$$        1,1,1;          ... white
% $$$        1,1,1];          % white
% $$$ % COMPUTE state occupancy
% $$$ spdur = cell2mat(cf(@(s,sts) sum(diff([s{sts,1}.data],1,2)),repmat({stc},[1,numel(statesThird)]),statesThird));
% $$$ % COMPUTE states' positions as blocks ordered left to right
% $$$ sppos = [0,cumsum(spdur(1:end-1))];
% $$$ 
% $$$ % PLOT Primary Level: all states mutually exclusive 
% $$$ for s = 1:numel(spdur);
% $$$     patch(sppos(s)+[0,spdur(s),spdur(s),0],... X
% $$$           [0,0,-1,-1],      ... Y
% $$$           sc(s,:));         % C
% $$$ end    
% $$$ 
% $$$ 
% $$$ 
% $$$ statesFourth = {'lloc&theta','lloc-theta','hloc&theta','hloc-theta',...
% $$$                'rear','lpause&theta','lpause-theta','hpause&theta','hpause-theta',...
% $$$                'groom','sit'};
% $$$ 
% $$$ sc  = [0,0,1;         ... l blue   lloc&theta
% $$$        0.5,0.5,1;     ... l violet lloc-theta
% $$$        0.25,0.25,0.75;... l cyan 
% $$$        0.5,0.5,1;     ... l violet lloc-theta
% $$$        1,1,1;         ... white    
% $$$        0,1,1;         ... cyan
% $$$        0,0.5,0.5;     ... dark cyan              
% $$$        0.25,0.75,0.75;... light cyan
% $$$        0,0.5,0.5;     ... dark cyan       
% $$$        1,1,1;         ... white
% $$$        1,1,1];          % white
% $$$ % COMPUTE state occupancy
% $$$ spdur = cell2mat(cf(@(s,sts) sum(diff([s{sts,1}.data],1,2)),repmat({stc},[1,numel(statesFourth)]),statesFourth));
% $$$ % COMPUTE states' positions as blocks ordered left to right
% $$$ sppos = [0,cumsum(spdur(1:end-1))];
% $$$ 
% $$$ % PLOT Primary Level: all states mutually exclusive 
% $$$ for s = 1:numel(spdur);
% $$$     patch(sppos(s)+[0,spdur(s),spdur(s),0],... X
% $$$           [-1,-1,-2,-2],      ... Y
% $$$           sc(s,:));         % C
% $$$ end    


set(gca,'YTick',[-1.5,-0.5,0.5,1.5]);
set(gca,'YTickLabels',{'HL theta','High/Low','theta','State'});
axis('tight');
set(gca,'XTick',[]);


% TESTING

configure_default_args();
MjgER2016_load_data();

[hfig,fig,fax,sax] = set_figure_layout(figure(666007),'A4','portrait',[],1.5,1.5,0.05,0.05);
hfig.PaperPositionMode = 'auto';
for t = 1:numel(Trials),
    subplot(numel(Trials),1,t);
    plot_state_durations(Trials{t});
    if   t ~= 1,set(gca,'YTickLabels',{});end
    hyl = ylabel(Trials{t}.filebase);
    hyl.Rotation = 0;
    hyl.FontSize = 8;
    hyl.HorizontalAlignment = 'right';
    hyl.VerticalAlignment = 'middle';
    hax = gca;
    hax.Position([1,3]) = [0.3,0.5];
end
hfig.CurrentAxes = fax;




FigName = ['stateDurations_',sessionListName,];
create_directory(fullfile(MTA_FIGURES_PATH,'behavior'));
print(gcf,'-depsc2',fullfile(MTA_FIGURES_PATH,'behavior',[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(MTA_FIGURES_PATH,'behavior',[FigName,'.png']));
