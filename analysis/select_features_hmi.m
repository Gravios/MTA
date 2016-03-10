function select_features_hmi(Trial,stc,fet,states)

[tstc,~,tfet] = resample_whole_state_bootstrap_noisy_trim(stc,fet,states);

fetInds = {};
gStates = states;

mixy = zeros([1+numel(gStates),tfet.size(2)]);
for s = 0:numel(gStates),
    
    if s==0,
        tStates = gStates;
    else
        tStates = gStates(~ismember(1:numel(gStates),s));
    end

    sm = stc2mat(tstc,tfet,tStates);
    [~,smat] = max(sm,[],2);
    smat(all(sm==0,2)) = 0;
    
    vind = smat&nniz(tfet);
    nind = sum(vind);

    for f = 1:tfet.size(2),

        if numel(tStates)==6,
            lsp{f} = mat2cell([prctile(tfet(vind,f),[1,99]),2^7],1,[1,1,1]);
        end
                
        edx = linspace(lsp{f}{:});
        edy = .5:numel(tStates)+.5;
        
        fx = tfet(vind,f);
        fy = smat(vind);

        [out,xb,yb,p]=hist2([fx,fy],edx,edy);
        pxy = out./nind;
        px = histc(fx,xb); px = px(1:end-1)/nind;
        py = histc(fy,yb); py = py(1:end-1)/nind;
        mixy(s+1,f) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
    end

end

hfig = figure(283823899);
for s = 1:numel(gStates),
    subplot(2,ceil(numel(gStates)/2),s);
    plot(mixy([1,s+1],:)');
    title(['MI of all states VS all without ' gStates{s}]);
    xlabel('Features');
    ylabel('Mutual Information (bits)');
end
reportfig(fullfile(Trial.path.data,'figures'),...
          hfig,...
          strjoin({Trial.filebase,fet.label,stc.mode,states{:}},'-'),...
          mfilename,...
          


figure,hold on,
dms = [];
for s = 1:numel(gStates)
    dm = (mixy(1,:)-mixy(s+1,:))';
    scatter(s,sum(dm(dm>0)),20);
    dms(end+1) = sum(dm(dm>0));
end

[~,sind] = max(dms);
dm = (mixy(1,:)-mixy(sind+1,:))';
%fetInds{end+1} =  find(dm>0.20);
fetInds{end+1} =  find(dm>0.20);

sind = 2; % rear
dm = (mixy(1,:)-mixy(sind+1,:))';
fetInds{end+1} =  find(dm>0.20);
