
sessionListName = 'MjgER2016';

FigDir = '/storage/gravio/figures/placefields';
mkdir(FigDir);


sessionList = get_session_list(sessionListName);
numTrials   = numel(sessionList);
sessionList = mat2cell(sessionList,1,ones([1,numTrials]));
overwrite = true;

for tind = 1:numTrials,
    
Trial = MTATrial.validate(sessionList{tind});
mkdir(fullfile(FigDir,Trial.filebase));

% LOAD behavioral state collection
stcMode = 'msnn_ppsvd_raux';
states  = {'loc','lloc','hloc','rear','pause','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
states{end+1} = 'theta';
statesCcg  = {'loc','lloc','hloc','rear','pause','lpause','hpause','theta'};
stc = label_bhv_reduced(Trial.load('stc','msnn_ppsvd'),Trial);
% load labeled behavior
try, Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.gnd','.stc.',stcMode,'.mat']);
catch err, disp(err);
    Trial.load('stc',[Trial.name,'.',Trial.maze.name,'.all','.stc.',stcMode,'.mat']);            
end

% LOAD
Trial.load('nq');
if numel(Trial.nq.eDist)~=size(Trial.spk.map,1),
    NeuronQuality(MTASession.validate(Trial.filebase),'overwrite',true);
    Trial.load('nq');    
end

 
pft = pfs_2d_theta(Trial,'overwrite',overwrite);
mrt = pft.maxRate;

% REDUCE clu list based on theta pfs max rate
units = select_units(Trial,18);
units = units(mrt(pft.data.clu(units))>1);
mrt = mrt(mrt>1);

% $$$ % COMPUTE place fields and subsampled estimate
% $$$ for sts = 1:numel(states),
% $$$     defargs = get_default_args('MjgER2016','MTAApfs','struct');
% $$$     defargs.units = units;
% $$$     defargs.states = states{sts};
% $$$     defargs.overwrite = overwrite;    
% $$$     defargs = struct2varargin(defargs);        
% $$$     pfbs{sts} = MTAApfs(Trial,defargs{:});      
% $$$ end

% COMPUTE place fields and subsampled estimate
for sts = 1:numel(states),
    defargs = get_default_args('MjgER2016','MTAAknnpfs','struct');
    defargs.units = units;
    defargs.states = states{sts};
    defargs.overwrite = overwrite;    
    defargs = struct2varargin(defargs);        
    pfk{sts} = MTAAknnpfs(Trial,defargs{:});      
end

% COMPUTE place fields and subsampled estimate
for sts = 1:numel(states),
    defargs = get_default_args('MjgER2016','MTAAknnpfs_bs','struct');
    defargs.units = units;
    defargs.states = states{sts};
    defargs.overwrite = overwrite;    
    defargs = struct2varargin(defargs);        
    pfk_bs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
end

% COMPUTE place fields and subsampled estimate
for sts = 1:numel(states),
    defargs = get_default_args('MjgER2016','MTAAknnpfs_perm','struct');
    defargs.units = units;
    defargs.states = states{sts};
    defargs.overwrite = overwrite;    
    defargs = struct2varargin(defargs);        
    pfk_perm{sts} = MTAAknnpfs_perm(Trial,defargs{:});      
end


% COMPUTE place fields and subsampled estimate
for sts = 1:numel(states),
    [bhvccg{sts},sper{sts}] = gen_bhv_ccg(Trial,statesCcg{sts},0.5);
end


xyz = Trial.load('xyz');
% Plot place fields 
% 1st row: Placefields 
% 2st row: Locations of onsets and offsets
% 3rd row: ccg onset of behavior 
% 4th row: ccg offset of behavior


hfig = figure();
hfig.Position = [423,         237,        1106,         447];
hfig.PaperPositionMode = 'auto';
u = 2;
for u = 1:numel(units),
    clf();
    unit = units(u);
    mpfsRate = max(cell2mat(cf(@(p,u) max(p.maxRate(u)),...
                               pfkbs,repmat({unit},[1,numel(bhvccg)]))));
    if mpfsRate<=0,mpfsRate=1;end % default to 1 if 0
    mccgRate = max(cell2mat(cf(@(c,u) max(max(c.ccg(:,u,:))),...
                               bhvccg,repmat({unit},[1,numel(bhvccg)]))));
    if mccgRate<=0,mccgRate=1;end % default to 1 if 0
    for s = 1:numel(states),

        % row 1 - PlaceFields
        subplot2(4,numel(states),1,s);
        plot(pfkbs{s},unit,[],[],mpfsRate);
        title(statesCcg{s});
        
        % row 2 - state positions
        subplot2(4,numel(states),2,s);
        plot(
        xlim([-500,500]);ylim([-500,500]);
        
        % row 2 - State Onset
        subplot2(4,numel(states),3,s);
        plot(bhvccg{s},unit,1);
        axis('tight');
        ylim([0,mccgRate]);
        title('onset');
        
        % row 3 - State Offset
        subplot2(4,numel(states),4,s);
        plot(bhvccg{s},unit,2);    
        axis('tight');    
        ylim([0,mccgRate]);
        title('offset');
        
    end
    suptitle(['Trial: ',Trial.filebase,' Unit: ',num2str(unit)]);

    FigName = ['pfs','_',Trial.filebase,'_unit-',num2str(unit)];
    %print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
    print(gcf,'-dpng',  fullfile(FigDir,Trial.filebase,[FigName,'.png']));
end


end
% $$$ 
% $$$ referenceTrial = 'Ed05-20140529.ont.all';
% $$$ pch = fet_HB_pitch(Trial);
% $$$ pch.map_to_reference_session(Trial,referenceTrial);    
% $$$ filter(pch,'ButFilter',3,1.5,'low');
% $$$ 
% $$$ figure();hold('on');
% $$$ eds = linspace(-pi/2,pi/2,50);
% $$$ ind = [Trial.stc{'lpause'}];
% $$$ hax = bar(eds,histc(pch(ind,3),eds),'histc');
% $$$ hax.FaceColor = 'b';
% $$$ hax.FaceAlpha = 0.5;
% $$$ hax.EdgeColor = 'b';
% $$$ hax.EdgeAlpha = 0.5;
% $$$ ind = [Trial.stc{'hpause'}];
% $$$ hax = bar(eds,histc(pch(ind,3),eds),'histc');
% $$$ hax.FaceColor = 'r';
% $$$ hax.FaceAlpha = 0.5;
% $$$ hax.EdgeColor = 'r';
% $$$ hax.EdgeAlpha = 0.5;


