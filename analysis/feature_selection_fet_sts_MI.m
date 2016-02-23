function  fsmi = feature_selection_fet_sts_MI(Trial,varargin);
[fet,stcMode,states,sampleRate] = DefaultArgs(varargin,{'fet_all','hand_labeled_rev2_jg',{'walk','rear','turn','pause','groom','sit'},12});

Trial = MTATrial.validate('jg05-20120317');
%Trial = MTATrial.validate(Trial);
stc = Trial.load('stc',stcMode);


RefTrial = MTATrial.validate('jg05-20120317');


fet = fet_all(Trial,sampleRate,RefTrial);
%fet = fet_tsne_rev15(Trial,sampleRate);
fet.data = [fet.data,fet.data.^2];
afet = fet.copy;
for sh = 1:fet.size(2);
    fet.data = [fet.data,circshift(afet.data',-sh)'.*afet.data];
end
fet = fet.unity;

[tstc,~,tfet] = resample_whole_state_bootstrap_noisy_trim(stc,fet,states);


fetE = MTADfet(Trial.spath,...
               [],...
               [],...
               sampleRate,...
               Trial.sync.copy,...
               Trial.sync.data(1),...
               [],'TimeSeries',[],'mi_Features','fet_mi_select','i');                  



fetInds = {};
stsOrd  = {};
gStates = states;

while numel(gStates) > 2,

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
fetInds{end+1} =  find(dm>0.50);



% $$$ fetS = fet.copy;
% $$$ fetS.data = fet(:,fetInds{end});
% $$$ aind = Trial.stc{['gper-' strjoin(stsOrd,'-')]}.cast('TimeSeries');
% $$$ aind.resample(fet);
% $$$ 
% $$$ [U,S,V] = svd(fetS(~~aind.data&nniz(fetS),:)'*fetS(~~aind.data&nniz(fetS),:));
% $$$ 
% $$$ svar = cumsum(diag(S)./sum(diag(S)));
% $$$ 
% $$$ fetE.data = [fetE.data,fetS.data*V(:,1:find(svar>.9,1,'first'))];
stsOrd{end+1} = gStates{sind};
gStates(sind) = [];

                 
end

% Used D to select final features
Dscore = abs(mean(tfet(tstc{gStates{1}},:))-mean(tfet(tstc{gStates{2}},:)))./sqrt(std(tfet(tstc{gStates{1}},:)).*std(tfet(tstc{gStates{2}},:)));

figure,plot(Dscore)

[~,dind] = sort(Dscore,'descend');
fetInds{end+1} = dind(1:20);


[smat] = stc2mat(Trial.stc,fetE,states);
net = patternnet(300);
%net.trainParam.showWindow = false;
%view(net);    
ind = nniz(fetE);
[net,tr] = train(net,fetE(ind,:)',~~smat(ind,:)');





%% Create Sit Feature
s = 6;

fetS = fet.copy;
fetS.data = fet(:,fetInds{1});
aind = Trial.stc{'a'}.cast('TimeSeries');
aind.resample(fet);

[U,S,V] = svd(fetS(~~aind.data&nniz(fetS),:)'*fetS(~~aind.data&nniz(fetS),:));
figure,plot(fetS.data*V(:,3));

px = fetS.copy;
fetE.data = fetS.data*V(:,1:3);
px.data = fetS.data;

eds = linspace(-20,40,100);
figure,hold on
ind = Trial.stc{'a-s'};
hs = bar(eds,histc(px(ind,1),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'s'};
hs = bar(eds,histc(px(ind,1),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';







%% Create Rear Feature
s = 2;
ind = find((mixy(1,:)-mixy(s+1,:))>.1);
[~,Rm,Rs] = nunity(fet(Trial.stc{'a-r'},ind));

fetR = fet.copy;
fetR.data = fet(:,ind);
fetR.unity([],Rm,Rs);
aind = Trial.stc{'a'};
[U,S,V] = svd(fetR(aind,:)'*fetR(aind,:));
pfigure,plot(fetR.data*V(:,1));

px = fetR.copy;
fetR.data = fetR.data*V(:,1);
px.data = fetR.data;



eds = linspace(-5,20,100);
figure,hold on
ind = Trial.stc{'a-r'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'r'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';


%% Create Sit Feature
s = 6;
ind = find((mixy(1,:)-mixy(s+1,:))>.2);
[~,Rm,Rs] = nunity(fet(Trial.stc{'a-s-r'},ind));

fetS = fet.copy;
fetS.data = fet(:,ind);
fetS.unity([],Rm,Rs);
aind = Trial.stc{'a-r'};
[U,S,V] = svd(fetS(aind,:)'*fetS(aind,:));
figure,plot(fetS.data*V(:,1));

px = fetS.copy;
fetS.data = fetS.data*V(:,1);
px.data = fetS.data;

eds = linspace(-20,30,100);
figure,hold on
ind = Trial.stc{'a-s'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'s'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';

%% Create groom Feature
s = 4;
ind = find((mixy(1,:)-mixy(s+1,:))>.1);
[~,Rm,Rs] = nunity(fet(Trial.stc{'a-s-r-m'},ind));

fetM = fet.copy;
fetM.data = fet(:,ind);
fetM.unity([],Rm,Rs);
aind = Trial.stc{'a-r-s'};
[U,S,V] = svd(fetM(aind,:)'*fetM(aind,:));
figure,plot(fetM.data*V(:,2));

fm = fetM.copy;


px = fetM.copy;
fetM.data = fm.data*V(:,1);
px.data = fetM.data;

eds = linspace(-20,20,100);
figure,hold on
ind = Trial.stc{'a-m'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'m'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';


%% Create groom Feature
s = 4;
ind = find((mixy(1,:)-mixy(s+1,:))>.1);
[~,Rm,Rs] = nunity(fet(Trial.stc{'a-s-r-m'},ind));

fetM = fet.copy;
fetM.data = fet(:,ind);
fetM.unity([],Rm,Rs);
aind = Trial.stc{'a-r-s'};
[U,S,V] = svd(fetM(aind,:)'*fetM(aind,:));
figure,plot(fetM.data*V(:,1));

fm = fetM.copy;


px = fetM.copy;
fetM.data = fm.data*V(:,1);
px.data = fetM.data;

eds = linspace(-20,20,100);
figure,hold on
ind = Trial.stc{'a-m-s-r'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'m'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';



%% Create Pause Feature
s = 3;
ind = find((mixy(1,:)-mixy(s+1,:))>.1);
[~,Rm,Rs] = nunity(fet(Trial.stc{'a-s-r-m-p'},ind));

fetP = fet.copy;
fetP.data = fet(:,ind);
fetP.unity([],Rm,Rs);
aind = Trial.stc{'a-r-s-m'};
[U,S,V] = svd(fetP(aind,:)'*fetP(aind,:));
figure,plot(fetP.data*V(:,));

fm = fetP.copy;


px = fetP.copy;
fetP.data = fm.data*V(:,6);
px.data = fetP.data;

eds = linspace(-20,20,100);
figure,hold on
ind = Trial.stc{'a-m-s-r'};
hs = bar(eds,histc(px(iond),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'c';
ind = Trial.stc{'m'};
hs = bar(eds,histc(px(ind),eds),'histc');
hs.FaceAlpha = .5;hs.FaceColor = 'r';





edxy = linspace(-15,30,100);
snd = Trial.stc{'a'};
figure,hist2([fetR(snd),fetS(snd)],edxy,edxy);
caxis([0,50])



figure,hold on,
for s = 1:numel(states)
    dm = (mixy(1,:)-mixy(s+1,:))';
    scatter(s,max(dm),20);
end



figure,hold on,
for s = 1:numel(states)
    dm = (mixy(1,:)-mixy(s+1,:))';
    scatter(s,mean(dm(dm>0))./(1+abs(mean(dm(dm<0)))),20);
    mean(dm(dm>0))./abs(mean(dm(dm<0)))
end

sm = stc2mat(stc,fet,states);
sm = [sum(sm,2)~=2,sum(sm,2)==2];
nind = any(~~sm,2);
net = patternnet(100);
net = train(net,fet(nind,dm>0)',sm(nind,:)');

figure,plot(min(bsxfun(@minus,mixy(2:end,:),mixy(1,:))'),'.')
