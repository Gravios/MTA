

train = true;
states = {'walk','rear','turn','pause','groom','sit'};
target = 'sit';
sampleRate = 12;



if train
    %Train Parm
    stcMode = 'hand_labeled_rev3_jg';
    Trial = MTATrial.validate('jg05-20120317');
    Trial.load('stc',stcMode);
    stc = Trial.stc.copy;
    RefTrial = [];
    rMean = []; rStd = [];
else
    %Test Parm
    %Trial = MTATrial.validate('Ed03-20140625');
    stcMode = 'hand_labeled_rev1_Ed';
    Trial = MTATrial.validate({'Ed05-20140529','all','ont'});    
    Trial.load('stc',stcMode);
    stc = Trial.stc.copy;
    RefTrial = MTATrial.validate('jg05-20120317');
    RefTrial.load('stc','hand_labeled_rev3_jg');
    rfet = fet_all(RefTrial,sampleRate,[]);
    rfet.data = [rfet.data,rfet.data.^2];
    rafet = rfet.copy;
    for sh = 1:rfet.size(2)-1;
        rfet.data = [rfet.data,circshift(rafet.data',-sh)'.*rafet.data];
    end
    [~,rMean,rStd] = unity(rfet);
    clear('rfet')
end


% LOAD all features
fet = fet_all(Trial,sampleRate,RefTrial);
fet.data = [fet.data,fet.data.^2];
afet = fet.copy;
for sh = 1:fet.size(2)-1;
    fet.data = [fet.data,circshift(afet.data',-sh)'.*afet.data];
end

% NORMALIZE features
if ~isempty(rMean)&&~isempty(rStd),
    fet = unity(fet,[],rMean,rStd);
else
    fet = fet.unity;
end


if train,
    [tstc,~,tfet] = resample_whole_state_bootstrap_noisy_trim(stc,fet,states);
    tstc.states{end}.data = [1,tfet.size(1)];
    [stateOrd,fetInds,miAstates] = select_features_hmi(Trial,tstc,tfet,states,false);

    mpn = struct;
    opn = struct;
    gStates = states;
% $$$     for sind = 1:numel(stateOrd),
% $$$         tstates = {[strjoin({gStates{find(cellfun(@isempty,regexp(gStates,['(',strjoin(stateOrd(1:sind),')|('),')'])))}},'+'),'&gper'],stateOrd{sind}};
% $$$         
% $$$ 
% $$$         % t-SNE 
% $$$         % sfet = fet.copy;
% $$$         % sfet.data = tfet(:,fetInds{sind});
% $$$         % mta_tsne(Trial,sfet,12,tstc,tstates,5,2,80,'ifReportFig',false,'overwrite',false);
% $$$ 
% $$$ 
% $$$         % NN labeling 
% $$$ % $$$         sfet = fet.copy;
% $$$ % $$$         sfet.data = fet(:,fetInds{sind});
% $$$ % $$$         sfet.label = [sfet.label '_' stateOrd{sind}]
% $$$ % $$$         [mpn.Stc,mpn.d_state,mpn.labelingStats, ...
% $$$ % $$$          mpn.labelingStatsMulti,mpn.model,mpn.p_state] = ...
% $$$ % $$$             bhv_nn_multi_patternnet(Trial,tstates,stc,sfet,[],[],200,3,'WSBNT');
% $$$        
% $$$         gStates(~cellfun(@isempty,regexp(gStates,['^',stateOrd{sind},'$'])))=[];
% $$$     end

afetinds = unique(vertcat(fetInds{:}));
sfet = fet.copy;
sfet.data = fet(:,afetinds);
sfet.label = [sfet.label '_MIselectedSubset'];
[mpn.Stc,mpn.d_state,mpn.labelingStats, ...
 mpn.labelingStatsMulti,mpn.model,mpn.p_state] = ...
    bhv_nn_multi_patternnet(Trial,states,stc,sfet,[],[],200,3,'WSBNT');


% Test if features without expansion work as well

oind = [repmat([1:59],1,2)',zeros([118,1])];
aind = oind(:,1);
for sh = 1:117,
    oind = [oind;[circshift(aind,-sh),aind]];
end


slind = oind(afetinds,:);
slind = oind(fetInds{6},:);
ofet = unique(reshape(slind,[],1));




% $$$ 
% $$$ for r = 1:4,
% $$$     %r = 1;
% $$$ hfig = figure;hold('on');
% $$$ eds = linspace(-0.25,0.75,50);
% $$$ tstates = states;
% $$$ tstates = {gStates{find(cellfun(@isempty,regexp(tstates,['(',strjoin(stateOrd(1:r-1),')|('),')'])))}};
% $$$ midsum = zeros([1,numel(tstates)]);
% $$$ for s = 1:(size(miAstates{r},1)-1),
% $$$     midiff = miAstates{r}(1,:)-miAstates{r}(1+s,:);
% $$$     midsum(s) = sum(abs(midiff));
% $$$     plot(eds,log10(histc(midiff,eds)+1),'linewidth',2);
% $$$ end
% $$$ xlabel(['States/feature MutInfo Difference between vs all vs all one ' ...
% $$$         'held out'])
% $$$ ylabel(['log10+1 count'])
% $$$ legend(cellfun(@strcat,tstates,repmat({' sum(abs(d(mi))):'},1,numel(tstates)),cellfun(@num2str,mat2cell(midsum,1,ones(1,numel(tstates))),'uniformoutput',false),'uniformoutput',false));
% $$$ 
% $$$ saveas(hfig,fullfile('/storage/gravio/manuscripts/man2015-jgEd-MoCap/Figures/Figure_3/',...
% $$$                 ['MIDIFF-' strjoin(tstates,'-') '.eps']),'epsc')

% $$$         reportfig(fullfile(Trial.path.data,'figures'),... 
% $$$                   hfig,...          Figure Handle 
% $$$                   strjoin({Trial.filebase,fet.label,stc.mode,tstates{:}},'-'),... FileName
% $$$                   'req20160310',...     Subdirectory 
% $$$                   false,...         Preview
% $$$                   strjoin({Trial.filebase,fet.label,stc.mode,tstates{:}},'-'),... Tag
% $$$                   strjoin({Trial.filebase,fet.label,stc.mode,tstates{:}},'-'),... Comment
% $$$                   200,...           Resolution
% $$$                   true,...          SaveFig
% $$$                   'png',...         Format
% $$$                   16,10);%          width & height (cm)

end
%'tsnem','mpn'
model = mpn.model;
model ='MTAC_BATCH-fet_all_MIselectedSubset_SR_12_NORM_0_REF_jg05-20120317.cof.all_STC_hand_labeled_rev3_jg_NN_200_NI_3_NN_multiPN_RAND_WSBNT';

sfet = fet.copy;
sfet.data = fet(:,afetinds);
sfet.label = [sfet.label '_MIselectedSubset'];
[opn.Stc,opn.d_state,opn.labelingStats, ...
 opn.labelingStatsMulti,opn.model,opn.p_state] = ...
    bhv_nn_multi_patternnet(Trial,states,stc,sfet,[],model,200,3,'WSBNT');
