
% Tasks 20180509
%
% 1. Compile: rate, size, position, and spatial information, of each place field state
% 2. plot: state x state patch distance 


% Compute shuffled state pfs
MjgER2016_load_data;



tind = 20; % jg05-20120312.cof.all
tind = 1; % jg05-20120312.cof.all

pfs = pfs_2d_states(Trials{tind},units{tind},[],[],false,'',false);
pfsArgs = struct('posShuffle',true);
pfsShuffled = pfs_2d_states(Trials{tind},units{tind},[],[],false,'',false,pfsArgs);

pft = pfs_2d_theta(Trials{tind},units{tind});
pftArgs = struct('posShuffle',true);
pftShuffled = pfs_2d_theta(Trials{tind},units{tind},false,true,pftArgs);


pfs = cat(2,{pft},pfs);
pfsShuffled = cat(2,{pftShuffled},pfsShuffled);

numStates = numel(pfs);

%pbins = interpParPfs.bins; ipp = interpParPfs;
pbins = pfs{1}.adata.bins; ipp = [];


bins = cell([1,2]);
[bins{:}] = ndgrid(pbins{:});

for u = 1:numel(units{tind});
clf();
mrate = max(cell2mat(cf(@(p,u) max(max(plot(p,u,'mean',false,[],true,0.5,false,[]))), pfs,repmat({units{tind}(u)},size(pfs)))));
for s = 1:numStates;

    rmapC = plot(pfs{s},units{tind}(u),'mean',false,[],true,0.5,false,ipp);
    rmapSm = plot(pfsShuffled{s},units{tind}(u),'mean',false,[],true,0.5,false,ipp);
    rmapSs = plot(pfsShuffled{s},units{tind}(u),'std',false,[],true,0.5,false,ipp);
% $$$     rmapSm = plot(pftShuffled,units{tind}(u),'mean',false,[],true,0.5,false,ipp);
% $$$     rmapSs = plot(pftShuffled,units{tind}(u),'std',false,[],true,0.5,false,ipp);
% $$$     rmapSm = plot(pfs{s},units{tind}(u),'mean',false,[],true,0.5,false,ipp);
% $$$     rmapSs = plot(pfs{s},units{tind}(u),'std',false,[],true,0.5,false,ipp);
    rmapSnr = plot(pfs{s},units{tind}(u),'snr',false,[],true,0.5,false,ipp);
    

    subplot2(numStates,5,s,1);  imagescnan({pbins{:},rmapC'},[0,mrate],[],true,[],[],[],@jet); 
    hold('on');                 contour(bins{:},rmapC,repmat(max(rmapC(:)),[1,2]).*0.5,'m')
                                title(num2str(units{tind}(u)));
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});
                            
    subplot2(numStates,5,s,2);  imagescnan({pbins{:},rmapSm'},[],[],true,[],[],[],@jet); 
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});    
    subplot2(numStates,5,s,3);  imagescnan({pbins{:},rmapSs'},[],[],true,[],[],[],@jet); 
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});


    zmap = (rmapC'-rmapSm')./rmapSs';
    subplot2(numStates,5,s,4);  imagescnan({pbins{:},zmap},[0,10],[],true,[],[],[],@jet);
    hold('on');                 contour(bins{:},zmap',repmat(3,[1,2]),'m')
                                title(pfs{s}.parameters.states);
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});
                               
    subplot2(numStates,5,s,5);  imagescnan({pbins{:},rmapSnr'},[0,10],[],true,[],[],[],@jet); 
    hold('on');                 contour(bins{:},rmapSnr,repmat(3,[1,2]),'m')
                                title(pfs{s}.parameters.states);
                                set(gca,'YTickLabels',{});set(gca,'XTickLabels',{});
end
waitforbuttonpress();
end




% $$$ cf(@(t,u,a) pfs_2d_states(t,u,[],[],false,'',false,a), Trials,units,repmat({pfsArgs},size(Trials)));
% $$$ cf(@(t,u,a) pfs_2d_theta(t,u,false,false,a),           Trials,units,repmat({pfsArgs},size(Trials)));

pfstats  = {};
pfbstats = {};
pfmstats = {};
for u = 1:numel(units{tind}),
    for s = 1:numStates,
        [pfstats{u,s},pfbstats{u,s},pfmstats{u,s}] = compute_placefield_stats(Trials{tind},...
                                                          pfs{s},units{tind}(u),3,true,pftShuffled);
    end    
end
