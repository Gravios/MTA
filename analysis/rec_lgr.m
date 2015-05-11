


%Sessions = SessionList('test_grp',...
%                '/storage/gravio/data/processed/xyz/',...
%                '/storage/gravio/data/processed/nlx/');

%% Model Info
train = true;
states = {'walk','rear','sit','turn','shake','groom'};
init_ns = numel(states);

Trial = MTATrial('jg05-20120317','all','cof');    
fet = fet_lgr(Trial);
fet.resample(30);
model_names = {};
State_cat_order = {};%repmat({''},[init_ns,1]);
lgrm = {};
pB = repmat({[]},[init_ns,1]);

 
for s = 1:init_ns,
    disp(['inter: ' num2str(s) ', finding best state'])
    temp_states =  states(cellfun(@isempty,regexp(states,['(',strjoin(State_cat_order,'|'),')'])));
    
    if ~isempty(State_cat_order),
        sws = ['-' strjoin(State_cat_order,'-')];
    else
        sws = '';
    end
    
    for i = 1:numel(temp_states),
        model_names(s,i) = {[Trial.filebase,'-','pop_lgr-' temp_states{i} sws]};%mfilename]};
        bhv_lgr(Trial,train,[temp_states(i),Trial.stc{['a-' temp_states{i} sws]}.label],...
                fet,model_names{s,i},false,false);
    end

    
    for i = 1:numel(temp_states),
        lgrm{s,i} = load([model_names{s,i} '-fet_lgr-model.mat']);
        pB{s} = cat(2,pB{s},lgrm{s,i}.B);
    end
    
    [~,BestStateInd] = max(max(abs(pB{s}(2:end,:))));
    State_cat_order{s} = temp_states{BestStateInd};
    
end


tstates = states;
tstates(2) =[];
d_prime = zeros([fet.size(2),numel(tstates)]);
for i = 1:numel(tstates)
    d_prime(:,i) = ((mean(fet(Trial.stc{tstates{i}},:))-mean(fet(Trial.stc{['gper-rear-' tstates{i}]},:)))...
        ./sqrt(.5*(var(fet(Trial.stc{tstates{i}},:))+var(fet(Trial.stc{['gper-rear-' tstates{i}]},:)))))';
end

wfet = MTADxyz('data',fet(Trial.stc{'w'},d_prime(:,1)>1),'sampleRate',fet.sampleRate);
awfet = MTADxyz('data',fet(Trial.stc{'a-w'},d_prime(:,1)>1),'sampleRate',fet.sampleRate);
figure,hist(wfet(randi(wfet.size(1),20000),1),10000);

ind = randi(wfet.size(1),[40000,1]);
figure,
subplot(211),bar(linspace(-.75,2,1000),histc(mean(wfet(ind,[1]),2),linspace(-.75,2,1000)),'histc');
subplot(212),bar(linspace(-.75,2,1000),histc(mean(awfet(ind,[1]),2),linspace(-.75,2,1000)),'histc');



%lgrm = cat(1,cell2mat(lgrm));
%figure
%figure,imagesc(lgrm(1).stats.coeffcorr'),colorbar,colormap jet

Trial = MTATrial('jg05-20120310');

Trial = MTATrial('jg05-20120317');
dsx = Trial.load('xyz');
vl = dsx.vel([2,7],[1,2]);
[ys,fs,ts] = fet_spec(Trial,vl,'mtchglong',false);
ys.resample(30);

ind = Trial.stc{'a'};
figure,
hist2([log10(clip(sum(ys(ind,1:30,1,1),2),.00001,300)),...
       log10(clip(sum(ys(ind,1:30,2,2),2),.00001,300))],...
       linspace(-3,2.5,100),...
       linspace(-3,2.5,100))

mfet = [log10(clip(sum(ys(ind,1:30,1,1),2),.00001,300)),...
        log10(clip(sum(ys(ind,1:30,2,2),2),.00001,300))];
   
idx = clusterdata(mfet,'linkage','ward','savememory','on','maxclust',4);

figure
hist2(mfet(idx==4,:),...
       linspace(-3,2.5,100),...
       linspace(-3,2.5,100))

vl.resample(ys);
vla = vl(ind,:);

aper = Trial.stc{'a',vl.sampleRate}.cast('TimeSeries');
aper.data = aper.data(1:end-1);
ida = nan(vl.size);
ida(logical(aper.data)) = idx;
figure,plot(log10(clip(vl.data,0.01,200)))
hold on,plot(ida-6)
               
               
               
                   
                   
      
   