function req20160310_5_genfigs(Trial)

Trial = MTATrial.validate(Trial);


dsd = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));
for s = 1:5,
    ds = load(fullfile(Trial.spath,['req20160310_4_accumStats',num2str(s),'.mat']));

    oind = [repmat([1:59],1,2)',zeros([118,1])];
    aind = oind(:,1);
    for sh = 1:117,
        oind = [oind;[circshift(aind,-sh),aind]];
    end
    slind = oind(dsd.fetInds{s},:);
    ofet =reshape(slind,[],1);
    best_inds = histc(ofet,1:59);
    [~,sbind] = sort(best_inds,'descend');
    [~,fetLInds] = sort([0;diff(cell2mat({ds.accum_acc}))],'descend');
    fetLInds(fetLInds==1)=[];
    fetLInds = [1;fetLInds];
    oFetInds = sbind(1:size(fetLInds));    
    bFetInds{s} = oFetInds(fetLInds);
end


save(fullfile(Trial.spath,['req20160310_5_genfigs.mat']),'bFetInds','-v7.3');


% $$$ 
% $$$ oind = [repmat([1:59],1,2)',zeros([118,1])];
% $$$ aind = oind(:,1);
% $$$ for sh = 1:117,
% $$$     oind = [oind;[circshift(aind,-sh),aind]];
% $$$ end
% $$$ 
% $$$ pobj = parpool('local');
% $$$ 
% $$$ accum_acc = zeros([round(afet.size(2)/2),numel(fetInds)]);
% $$$ accum_pre = zeros([round(afet.size(2)/2),numel(fetInds)]);
% $$$ accum_sen = zeros([round(afet.size(2)/2),numel(fetInds)]);
% $$$ for s = 1:numel(fetInds);
% $$$     slind = oind(fetInds{s},:);
% $$$     ofet =reshape(slind,[],1);
% $$$     best_inds = histc(ofet,1:59);
% $$$     [~,sbind] = sort(best_inds,'descend');
% $$$ 
% $$$     for f = 1:round(numel(sbind)/2),    
% $$$         opn = struct;
% $$$         sub_fet = afet.copy;
% $$$         sub_fet.data = afet(:,sbind(1:f));
% $$$         sub_fet.label = [afet.label '-req20160310-' stateOrd{1} '-' num2str(f)];
% $$$         sub_fet.key = 'x';
% $$$         sub_fet.updateFilename(Trial);
% $$$         
% $$$         model = ['MTAC_BATCH-' stateOrd{s} sub_fet.label ...
% $$$                  '_SR_'  num2str(sub_fet.sampleRate) ...
% $$$                  '_NORM_' num2str(0) ...             
% $$$                  '_REF_' Trial.filebase ...
% $$$                  '_STC_' stcMode ...
% $$$                  '_NN_'  num2str(nNeurons)...
% $$$                  '_NI_'  num2str(nIter)...
% $$$                  '_'     'NN_multiPN'...
% $$$                  '_'     'RAND_' rndMethod];
% $$$ 
% $$$         [opn.Stc,opn.d_state,opn.labelingStats, ...
% $$$          opn.labelingStatsMulti,opn.model,opn.p_state] = ...
% $$$             bhv_nn_multi_patternnet(Trial,states,stc,sub_fet,...
% $$$                                     [], model,nNeurons,nIter,...
% $$$                                     rndMethod,'targetState',stateOrd{s});
% $$$ 
% $$$         accum_acc(f,s) = opn.labelingStats.accuracy;
% $$$         accum_pre(f,s) = opn.labelingStats.precision(1);
% $$$         accum_sen(f,s) = opn.labelingStats.sensitivity(2);
% $$$     end
% $$$     
% $$$ end
% $$$ 
% $$$ delete(pobj)
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


