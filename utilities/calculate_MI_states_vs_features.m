function calculate_MI_states_vs_features(stc,fet,states)
% function mixy = calculate_MI_states_vs_features(stc,fet,states)
% MI:=Mutual Information
%
% The first row contains the MI between state labels and features.
% The subsequent rows contains MI between state labels with one
% held out and features
% 

mixy = zeros([1+numel(states),fet.size(2)]);
for s = 0:numel(states),
    
    if s==0,
        tStates = states;
    else
        tStates = states(~ismember(1:numel(states),s));
    end

    sm = stc2mat(stc,fet,tStates);
    [~,smat] = max(sm,[],2);
    smat(all(sm==0,2)) = 0;
    
    vind = smat&nniz(fet);
    nind = sum(vind);

    for f = 1:fet.size(2),

        if numel(tStates)==6,
            lsp{f} = mat2cell([prctile(fet(vind,f),[1,99]),2^7],1,[1,1,1]);
        end
                
        edx = linspace(lsp{f}{:});
        edy = .5:numel(tStates)+.5;
        
        fx = fet(vind,f);
        fy = smat(vind);

        [out,xb,yb,p]=hist2([fx,fy],edx,edy);
        pxy = out./nind;
        px = histc(fx,xb); px = px(1:end-1)/nind;
        py = histc(fy,yb); py = py(1:end-1)/nind;
        mixy(s+1,f) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
    end

end
