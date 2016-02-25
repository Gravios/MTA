

states = {'walk','rear','turn','pause','groom','sit'};
sampleRate = 12;

%Train Parm
Trial = MTATrial.validate('jg05-20120317');
Trial.load('stc','hand_labeled_rev3_jg');
RefTrial = [];
rMean = []; rStd = [];

%Test Parm
Trial = MTATrial.validate('Ed03-20140625');
stc = Trial.load('stc','hand_labeled_rev1_Ed');
RefTrial = MTATrial.validate('jg05-20120317');
RefTrial.load('stc','hand_labeled_rev3_jg');
[~,rMean,rStd] = unity(fet_all(RefTrial,sampleRate,[]));

fet = fet_all(Trial,sampleRate,RefTrial);
fet.data = [fet.data,fet.data.^2];
afet = fet.copy;
for sh = 1:fet.size(2)-1;
    fet.data = [fet.data,circshift(afet.data',-sh)'.*afet.data];
end

if ~isempty(rMean)&&~isempty(rStd),
    fet = unity(fet,[],rMean,rStd);
else
    fet = fet.unity;
end


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

figure,
for s = 1:numel(gStates)
    subplot(2,3,s);
    plot(mixy([1,s+1],:)');
end

sind = 2; % rear
dm = (mixy(1,:)-mixy(sind+1,:))';
fetInds{end+1} =  find(dm>0.20);


smat = [[ones(7500,1);zeros(7500,1);ones(30000,1)],[zeros(7500,1);ones(7500,1);zeros(30000,1)]];
[net,tr] = train(net,tfet(nniz(tfet),:)',smat(nniz(tfet),:)');

d_state = net(fet.data')';

stc = Trial.load('stc','hand_labeled_rev3_jg');
