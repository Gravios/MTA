function  fsmi = feature_selection_fet_sts_MI(Trial,varargin);
[fet,stcMode,states,sampleRate] = DefaultArgs(varargin,{'fet_all','hand_labeled_rev2_jg',{'walk','rear','turn','pause','groom','sit'},12});

Trial = MTATrial.validate(Trial);

stc = Trial.load('stc',stcMode);

fet = fet_all(Trial,sampleRate);
[tstc,~,tfet] = resample_whole_state_bootstrap(stc,fet,states);



%states = {'walk','rear','turn','pause','groom','sit'};
%states = {'walk','turn','pause','groom','sit'};
%states = {'walk','turn','pause','groom'};
%states = {'rear','turn','sit'};

mixy = zeros([1+numel(states),fet.size(2)]);
for s = 0:numel(states),
    
    if s==0,
        tStates = states;
    else
        tStates = states(~ismember(1:numel(states),s));
    end

    sm = stc2mat(tstc,tfet,tStates);
    [~,smat] = max(sm,[],2);
    smat(all(sm==0,2)) = 0;
    %vind = Trial.stc{'a'}.cast('TimeSeries');
    %vind.resample(fet);
    %vind = logical(vind.data)&
    vind = smat&nniz(tfet);
    nind = sum(vind);

    for f = 1:tfet.size(2),

        lsp = mat2cell([prctile(tfet(vind,f),[1,99]),2^7],1,[1,1,1]);
        
        edx = linspace(lsp{:});
        edy = 1:numel(tStates);
        
        fx = tfet(vind,f);
        fy = smat(vind);

        [out,xb,yb,p]=hist2([fx,fy],edx,edy);
        pxy = out./nind;
        px = histc(fx,xb); px = px(1:end-1)/nind;
        py = histc(fy,yb); py = py(1:end-1)/nind;
        mixy(s+1,f) = nansum(nansum(pxy.*log2(pxy./(px*py'))));
    end

end

figure,
for s = 1:numel(states)
    subplot(2,3,s);
    plot(mixy([1,s+1],:)');
end

figure,hold on,
for s = 1:numel(states)
    dm = (mixy(1,:)-mixy(s+1,:))';
    scatter(s,sum(dm(dm>0)),20);
end



figure,hold on,
for s = 1:numel(states)
    dm = (mixy(1,:)-mixy(s+1,:))';
    scatter(s,max(dm),20);
end



figure,hold on,
for s = 1:numel(states)
    dm = (mixy(1,:)-mixy(s+1,:))';
    scatter(s,mean(dm(dm>0))./abs(mean(dm(dm<0))),20);
    mean(dm(dm>0))./abs(mean(dm(dm<0)))
end

sm = stc2mat(stc,fet,states);
sm = [sum(sm,2)~=2,sum(sm,2)==2];
nind = any(~~sm,2);
net = patternnet(100);
net = train(net,fet(nind,dm>0)',sm(nind,:)');

figure,plot(min(bsxfun(@minus,mixy(2:end,:),mixy(1,:))'),'.')
