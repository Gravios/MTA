function req20151007_A1(Trial,varargin)
[figPath,sampleRate,sts,stc_mode,initial_dims,perplexity] = ...
    DefaultArgs(varargin,{'/storage/gravio/figures/',15,{'walk','rear'},'auto_wbhr',5,100});

SAVEFIG = true;
osts = numel(sts)+1;


%% Load the State Collection
try
    Trial.stc.load(Trial,stc_mode);
catch
    if strcmp(stc_mode,'auto_wbhr'),        
        Trial = labelBhv(Trial);
        Trial.stc.load(Trial,'auto_wbhr');
    else
        error('MTAStateCollection:ModeDoesNotExist');
    end
end

%Check if desired states exist within the state collection (Trial.stc)
try 
    for s = 1:numel(sts),
        Trial.stc{sts{s}};
    end
catch
    error('MTAStateCollection:StateDoesNotExist');
end


%%Load features
[fet,fett,fetd] = fet_tsne(Trial,sampleRate,true);


[asmat,labels,keys] =  stc2mat(Trial.stc,fet,sts);
asmat = cat(2,asmat,~any(asmat,2));
asmat = MTADxyz('data',asmat,'sampleRate',fet.sampleRate);
[~,asmat.data] = max(asmat.data(:,[1:osts-1,end]),[],2);
c = jet(numel(unique(asmat)));
csmat = asmat.copy; 
csmat.data = c(csmat.data,:);

ind = Trial.stc{'a'};

mfet = nunity(fet(ind,:));
msmat = csmat(ind,:);

start = 1;
skip = 2;
stop = size(mfet,1);
no_dims = 2;

ind = start:skip:stop;

figTitle = ['tSNE-sampleRate_' num2str(sampleRate) '-ind_' num2str(start) '_' ...
            num2str(skip) '_' num2str(stop) '-perplexity_' ...
            num2str(perplexity) '-initial_dims_' num2str(initial_dims) ...
            '-no_dims_' num2str(no_dims)];

fileTSNE = fullfile(Trial.spath,[Trial.filebase '-req20150925-' figTitle '.mat']);

if  exist(fileTSNE,'file');
    ds = load(fileTSNE);
    mappedX = ds.mappedX;
else
    mappedX = tsne(mfet(ind,:), msmat(ind,:), no_dims, initial_dims, perplexity);
    save(fileTSNE,'mappedX');
end
    


hfig = figure(3923923);clf,hold on
mc = msmat(ind,:);
for nc = 1:osts,
    nind = all(bsxfun(@eq,c(nc,:),mc),2);
    h = scatter(mappedX(nind,1),mappedX(nind,2),10,mc(nind,:));
    h.MarkerFaceColor = h.CData(1,:);
end
legend([sts(1:osts-1),{'unlabeled'}])
xlim([min(mappedX(:,1))-5,max(mappedX(:,1))+30]);
ylim([min(mappedX(:,2))-5,max(mappedX(:,2))+5]);
daspect([1,1,1])


reportfig(figPath, hfig, 'tsne', 'req', false,Trial.filebase,figTitle,[],SAVEFIG,'png',[],[],[],true);


mtfet =  MTADxyz('data',mfet(ind,:),'sampleRate',fet.sampleRate);
mtpos =  MTADxyz('data',mappedX,'sampleRate',fet.sampleRate);

[RateMap,Bins,distdw,distIndw]= PlotKNNPF(Trial,mtfet,mtpos,[5,5],20,5,'xy',[],[],[-110,110;-110,110]);
rmap = reshape(RateMap,numel(Bins{1}),numel(Bins{2}),[]);



hfig = figure(38381);
for i = 1:size(rmap,3);
    clf
    [ha,hc] = imagescnan({Bins{1},Bins{2},rmap(:,:,i)},[0,3],false, ...
                        true,[0,0,0]);
    ylabel(hc,'mean z-score','FontName','Courier','FontSize',12);
    axis xy
    title(fett{i},'FontName','Courier','FontSize',12)
    daspect([1,1,1])
   
    reportfig(figPath, hfig, 'tsne', 'req',false,fett{i},fetd{i},[],SAVEFIG,'png',8,8);
end
