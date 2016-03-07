
train = false;
test = ~train;
states = {'walk','rear','turn','pause','groom','sit'};
target = 'sit';
sampleRate = 12;
gStates = states;


if train
    %Train Parm
    Trial = MTATrial.validate('jg05-20120317');
    Trial.load('stc','hand_labeled_rev3_jg');
    stc = Trial.stc.copy;
    RefTrial = [];
    rMean = []; rStd = [];
end

if test
    %Test Parm
    %Trial = MTATrial.validate('Ed03-20140625');
    Trial = MTATrial.validate({'Ed05-20140529','all','ont'});    
    Trial.load('stc','hand_labeled_rev1_Ed');
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


if train
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

    %state matrix of union between all states minus trarget and target
    smat = stc2mat(tstc,tfet,{[strjoin({states{find(cellfun(@isempty,regexp(states,states{sind})))}},'+'),'&gper'],states{sind}});
    net = patternnet(200);
    [net,tr] = train(net,tfet(nniz(tfet),fetInds{end})',smat(nniz(tfet),:)');

end



asmat = MTADfet('data',tstc,'sampleRate',fet.sampleRate);
[~,asmat.data] = max(asmat.data,[],2);
c = jet(numel(states));
c = [0,0,1;...
     1,0,0;...
     0,1,0;...
     0,1,1;...
     1,0,1;...
     1,1,0;];
csmat = asmat.copy; 
csmat.data = c(csmat.data,:);


osts = numel(states);
hfig = figure(3923924);clf
hold on;
mc = csmat(ind,:);
for nc = 1:osts,
    nind = all(bsxfun(@eq,c(nc,:),mc),2);
    h = scatter(mappedX(nind,1),mappedX(nind,2),2,mc(nind,:));
    try,h.MarkerFaceColor = h.CData(1,:);end
end
legend(states,'location','south','Orientation','horizontal');




d_state = net(fet(:,fetInds{end})')';

figure,plot(d_state)
Lines(stc{'r',12}(:),[],'r');

rper = stc{'r',12};
rper.cast('TimeSeries');
rper.resample(fet);
rper.data(isnan(rper.data)) = 0;

aper = stc{'a',12};
aper.cast('TimeSeries');
aper.resample(fet);
aper.data(isnan(aper.data)) = 0;

dper = stc{'a',12};
dper.data = ThreshCross(d_state(:,2),0.5,3);
dper = dper+[-.1,.2];
dper.cast('TimeSeries');
dper.resample(fet);
dper.data(isnan(dper.data)) = 0;

sum(aper.data.*(rper.data&dper.data))./sum(aper.data.*rper.data)
sum(aper.data.*(rper.data&dper.data))./sum(aper.data.*dper.data)

figure
hold on,plot(d_state(:,2)) 
Lines(stc{'r',12}(:),[],'r');