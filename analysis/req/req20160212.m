





NEW_SAMPLE_RATE = 30;
HOST_PATH = '/storage/gravio/figures/'; % Where reportfig should save stuff
states = {'walk','rear','turn','pause','groom','sit'};
featureSet = 'fet_tsne_rev15';
normalize = true;
map2reference = true;
sessionSet = 'hand_labeled_Ed';
mfilename = 'req20151203';


%Reference Trial Stuff
RefTrial = MTATrial('jg05-20120317');
if normalize,
    RefState = RefTrial.stc{'a'};
    rfet = feval(featureSet,RefTrial,NEW_SAMPLE_RATE);
    rfet.data = rfet(RefState,:);
    [~,Rmean,Rstd] = unity(rfet,[],[],[],[]);
end


slist = get_session_list(sessionSet);

sfet = [];
Stc = {};
sts = [];
for s = 1:numel(slist),
    Trial = MTATrial(slist(s).sessionName,slist(s).trialName,slist(s).mazeName);
    Trial.load('stc',slist(s).stcMode);
    Stc = Trial.stc.copy;
    % Load Feature matrix of the session    
    [tfet,fett,fetd] = feval(featureSet,Trial,NEW_SAMPLE_RATE);

    if strcmp(Trial.filebase,RefTrial.filebase)&&map2reference,
        tfet.map_to_reference_session(Trial,RefTrial);
    end

    if s == 1,
        fet = tfet.copy;
        fet.data = fet(Stc{'a'},:);
    else
        fet.data = cat(1,fet.data,tfet(Stc{'a'},:));
    end
    tsts = [];
    for t = states,
        tper = Stc{t{1}};
        tper.cast('TimeSeries');
        tper.resample(fet.sampleRate);
        tsts = cat(2,tsts,tper(Stc{'a'}));
    end     
    sts = cat(1,sts,tsts);
    sts(~nniz(fet),:) = [];        
    fet.data(~nniz(fet),:) = [];
    sfet(end+1) = length(fet.data);


end


edx = linspace(40,110,100);
figure,
bar(edx,histc(fet(:,1),edx),'histc')

edy = linspace(56,120,100);
figure,
bar(edy,histc(fet(:,2),edy),'histc')

edy = linspace(-3,2,100);
figure,
bar(edy,histc(fet(:,5),edy),'histc')

edy = linspace(0,1.5,100);
figure,
bar(edy,histc(fet(:,8),edy),'histc')

edy = linspace(-1,1.5,100);
figure,
bar(edy,histc(fet(:,9),edy),'histc')


[msh,x] = MakeUniformDistr(fet(:,1),40,110);edx = linspace(40,110,100);
[msu,x] = MakeUniformDistr(fet(:,2),56,120);
edx = linspace(56,120,100);

[msv,x] = MakeUniformDistr(fet(:,5),-3,2);  
edy = linspace(-3,2,100); 

[msb,x] = MakeUniformDistr(fet(:,6),-3,2);  edy = linspace(-3,2,100); 
[msa,x] = MakeUniformDistr(fet(:,8),0,1.5); edy = linspace(0,1.5,100);
[msa,x] = MakeUniformDistr(fet(:,9),-0.8,1.5); edy = linspace(-0.8,1.5,100);

figure,
for s = 1:6
subplot(2,3,s)
ind = ~~sts(:,s);
hist2([msu(ind),msv(ind)],edx,edy)
caxis([0,50])
end