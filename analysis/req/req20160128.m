function [out] = req20160128(varargin)
%
% DEFARGS ----------------------------------------------------------------------
%defargs = struct('fetSet',   'fet_bref',                                     ...
%                 'mode',     'display',                                      ...
%                 'tag_preprocessing',   '+seh+',                             ...
%                 'tag_postprocessing',  '',                                  ...
%                 'rlist',    'training_hand_labeled',                        ...
%                 'slist',     {{'hand_labeled_jg';'hand_labeled_Ed'}},       ...
%                 'sampleRate',15,                                            ...
%                 'nNeurons',  200,                                           ...
%                 'nIter',     100,                                           ...
%                 'states',    {{'walk','rear','turn','pause','groom','sit'}},...
%                 'rndMethod', 'WSBNT',                                       ...
%                 'norm',      true,                                          ...
%                 'mref',      true,                                          ...
%                 'prctTrain', []                                             ...
% );%-----------------------------------------------------------------------------

% DEFARGS ----------------------------------------------------------------------
defargs = struct('fetSet',   'fet_bref',                                     ...
                 'mode',     'display',                                      ...
                 'tag_preprocessing',   '+seh+',                             ...
                 'tag_postprocessing',  '',                                  ...
                 'rlist',    'training_hand_labeled',                        ...
                 'slist',     {{'hand_labeled_jg';'hand_labeled_Ed'}},       ...
                 'sampleRate',15,                                            ...
                 'nNeurons',  200,                                           ...
                 'nIter',     100,                                           ...
                 'states',    {{'walk','rear','turn','pause','groom','sit'}},...
                 'rndMethod', 'WSBNT',                                       ...
                 'norm',      true,                                          ...
                 'mref',      true,                                          ...
                 'prctTrain', [],                                            ...
                 'sessionList','hand_labeled'                                ...
);%-----------------------------------------------------------------------------


[fetSet, mode, tag_preprocessing, tag_postprocessing,rlist, slist, sampleRate,...
 nNeurons, nIter, states, rndMethod, norm, mref,prctTrain,sessionList] = ...
    DefaultArgs(varargin,defargs,'--struct');

rlist = get_session_list(rlist);

out = [];


% MAIN -------------------------------------------------------------------------

switch mode
  case 'train' %Trian bhv_nn_multi_patternnet 
    for s = rlist

        clear('mod');
        mod.states     = states;
        mod.stcMode    = s.stcMode;
        mod.featureSet = fetSet;
        mod.sampleRate = sampleRate;
        mod.nNeurons   = nNeurons;
        mod.nIter      = nIter;
        mod.randomizationMethod = rndMethod;
        mod.normalize = norm;
        mod.map2reference = mref;        
        mod.tag = tag_preprocessing;
        mod.prctTrain = prctTrain;
        mod.sessionList = sessionList;

        argin = struct2varargin(mod);
        
        bhv_nn_multi_patternnet(MTATrial.validate(s),argin{:});
    end



  case 'compute' % labels for trials in slist
    % Ed 
    for sli = 1:numel(slist),
        for rli = 1:numel(rlist),
            SesList = get_session_list(slist{sli});

            if ~isempty(prctTrain), prctTrainTag = ['_PRT_',prctTrain]; else, prctTrainTag = ''; end

            model = ['MTAC_BATCH-' tag_preprocessing fetSet ...
                     '_SR_' num2str(sampleRate)             ...
                     '_NORM_' num2str(norm)                 ...
                     '_REF_' rlist(rli).sessionName,        ...
                             '.' rlist(rli).mazeName        ...
                             '.' rlist(rli).trialName       ...
                     '_STC_' rlist(rli).stcMode             ...
                     '_NN_' num2str(nNeurons)               ...
                     '_NI_' num2str(nIter)                  ...
                     prctTrainTag                           ...                     
                     '_NN_multiPN_RAND_' rndMethod];

            stc = {}; d_state = {};p_state = {}; ls = {}; lsm = {};mdl = {};
            for s = SesList
                Trial = MTATrial.validate(s);
                Trial.load('stc',s.stcMode);
                clear('mod');
                mod.states     = states;
                mod.stcMode    = s.stcMode;
                mod.featureSet = fetSet;
                mod.model      = model;
                mod.sampleRate = sampleRate;
                mod.nNeurons   = nNeurons;
                mod.nIter      = nIter;
                mod.map2reference = mref;
                mod.normalize = norm;
                mod.prctTrain = prctTrain;
                mod.sessionList = sessionList;
                argin = struct2varargin(mod);
                [stc{end+1},d_state{end+1},ls{end+1},lsm{end+1},mdl{end+1},p_state{end+1}] = bhv_nn_multi_patternnet(Trial,argin{:});
                stc{end}.save(1);
            end

            if mod.map2reference,
                mapped = '-map2ref';
            else
                mapped = '';
            end
            
            save(fullfile(MTASession().path.project,'analysis',[slist{sli},'-',model,mapped,'.mat']),...
                 '-v7.3','slist','rlist','nNeurons','nIter','sampleRate','model','fetSet','rndMethod',...
                 'states','stc','d_state','p_state','ls','mdl');

        end
    end




  case 'display'
    if ~isempty(prctTrain), prctTrainTag = ['_PRT_',prctTrain]; else, prctTrainTag = ''; end                        
    for sli = 1:numel(slist),
        for rli = 1:numel(rlist),    


            SesList = get_session_list(slist{sli});
            SesList = {SesList(:).sessionName};
            
            nses = numel(SesList);
            nsts = numel(states);


            model = ['MTAC_BATCH-' tag_preprocessing fetSet ...
                     '_SR_' num2str(sampleRate) ...
                     '_NORM_' num2str(norm) ...         
                     '_REF_' rlist(rli).sessionName, '.' rlist(rli).mazeName '.' rlist(rli).trialName ...
                     '_STC_' rlist(rli).stcMode ...
                     '_NN_' num2str(nNeurons) ...
                     '_NI_' num2str(nIter) ...   
                     prctTrainTag                           ...                                          
                     '_NN_multiPN_RAND_' rndMethod];


            if mref
                mapped = '-map2ref';
            else
                mapped = '';
            end
            
            
            load(fullfile(MTASession().path.project,'analysis',[slist{sli},'-',model,...
                                tag_postprocessing,mapped,'.mat']))
            %slist = {'hand_labeled_jg';'hand_labeled_Ed'};


            %figure
            hfig = figure(3923992);clf
            set(hfig,'Position',       [40         40        1500         500])

            prop = 'accuracy';

            subplot(131);
            plot(cell2mat(cellfun(@subsref,ls, ...
                                  repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))'.*100,'d')
            xlim([0,5])
            ylim([0,100])
            ylabel(prop);
            title({['Training Set: ' rlist(rli).sessionName],...
                   ['Randomization: ' rndMethod ],...
                   ['Labeling Set: ' slist{sli}],...
                   ['Feature  Set: ' fetSet],...
                   ['postProc: ' tag_postprocessing]});
            set(gca,'XTickLabelMode','manual');
            set(gca,'XTick',1:numel(SesList));
            set(gca,'XTickLabel',SesList);
            set(gca,'XTickLabelRotation',90);
            pause(.1)

            prop = 'precision';
            subplot(132);plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                                                       repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',nses,nsts)','d-')
            xlim([0,7])
            hax = gca;
            hax.XTickLabelMode = 'manual';
            hax.XTickLabel = cat(2,{''},states,{''});
            ylim([0,100])
            ylabel(prop)
            title({['Training Set: ' rlist(rli).sessionName],...
                   ['Randomization: ' rndMethod ],...
                   ['Labeling Set: ' slist{sli}],...
                   ['Feature  Set: ' fetSet],...
                   ['postProc: ' tag_postprocessing]});
            legend(SesList,'location','southwest');
            pause(.1)

            prop = 'sensitivity';

            subplot(133);plot(reshape(cell2mat(cellfun(@subsref,ls, ...
                                                       repmat({substruct('.',prop)},[1,numel(ls)]),'uniformoutput',false))',nsts,nses),'d-')
            xlim([0,7])
            hax = gca;
            hax.XTickLabelMode = 'manual';
            hax.XTickLabel = cat(2,{''},states,{''});
            ylim([0,100])
            ylabel(prop)
            title({['Training Set: ' rlist(rli).sessionName],...
                   ['Randomization: ' rndMethod ],...
                   ['Labeling Set: ' slist{sli}],...
                   ['Feature  Set: ' fetSet],...
                   ['postProc: ' tag_postprocessing]});
            legend(SesList,'location','southwest');
            pause(3)



% $$$             reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$             hfig,                                ... Figure handle
% $$$             [mfilename],                         ... Figure Set Name
% $$$             'req',                               ... Directory where figures reside
% $$$             false,                               ... Do Not Preview
% $$$             ['Fet:',fetSet,'Ref:',rlist(rli).sessionName,'.',...
% $$$                                   rlist(rli).mazeName,'.',...
% $$$                                   rlist(rli).trialName,...
% $$$              ' L:',slist{sli}],...thmb_cap
% $$$             ['Fet:',fetSet,...
% $$$              'Ref:',rlist(rli).sessionName,'.',...
% $$$                     rlist(rli).mazeName,'.',...
% $$$                     rlist(rli).trialName,...
% $$$              ' L:',slist{sli}],...exp_cap
% $$$             [],                                  ... Resolution
% $$$             false,                               ... Do Not Save FIG
% $$$             'png',9,4);                                % Output Format
            
        end
    end
    
    
end



% $$$ 
% $$$ figure,hold on,
% $$$ fetSet  = 'fet_tsne_rev15';
% $$$ k = 1;
% $$$ htl = {};
% $$$ set(gcf,'interpreter','none')
% $$$ for sli = 1:numel(slist),
% $$$     for rli = 1:numel(refTrial),    
% $$$         prop = 'accuracy';
% $$$ 
% $$$         model = ['MTAC_BATCH-' fetSet ...
% $$$                  '_SR_' num2str(sampleRate) ...
% $$$                  '_NORM_' num2str(norm) ...
% $$$                  '_REF_' rlist(rli).sessionName, '.' ...
% $$$                  rlist(rli).mazeName '.' ...
% $$$                  rlist(rli).trialName ...
% $$$                  '_STC_' rlist(rli).stcMode ...
% $$$                  '_NN_' num2str(nNeurons) ...
% $$$                  '_NI_' num2str(nIter) ...
% $$$                  '_NN_multiPN_RAND_' rndMethod];
% $$$ 
% $$$         ds = load(fullfile(MTASession().path.project,'analysis',...
% $$$                            [slist{sli},'-',model,mapped,'.mat']));
% $$$ 
% $$$         scat = scatter(k*ones([numel(ds.ls),1]),cell2mat(cellfun(@subsref,ds.ls, ...
% $$$                                                           repmat({substruct('.',prop)},[1,numel(ds.ls)]),'uniformoutput',false))'.*100,30);
% $$$ 
% $$$         scat.MarkerFaceColor = scat.CData;
% $$$         k=k+1;
% $$$         htl = cat(2,htl,{['{' slist{sli} ' - ' rlist(rli).sessionName '}']});
% $$$     end
% $$$ end
% $$$ 
% $$$ xlim([0.5,k-0.5]);
% $$$ ylim([50,100])
% $$$ 
% $$$ set(gca,'XTickLabelMode','manual');
% $$$ set(gca,'XTick',1:k-1);
% $$$ set(gca,'XTickLabel',{});
% $$$ 
% $$$ 
% $$$ set(gca,'XTickLabel',htl);
% $$$ set(gca,'XTickLabelRotation',90);
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
