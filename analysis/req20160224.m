
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
    hmi = select_features_hmi(Trial,stc,fet,states)
    %state matrix of union between all states minus trarget and target

    smat = stc2mat(tstc,tfet,{[strjoin({states{find(cellfun(@isempty,regexp(states,states{sind})))}},'+'),'&gper'],states{sind}});

    
    net = patternnet(200);
    [net,tr] = train(net,tfet(nniz(tfet),fetInds{end})',smat(nniz(tfet),:)');

    mta_tsne(Trial,tfet,
end





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