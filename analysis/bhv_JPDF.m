function bhv_JPDF(Trial,var1,var2,varargin)

ind = Trial.stc{'a'};

[v1nbins,v2nbins,v1lim,v2lim,v1Label,v2Label,states,tag] = DefaultArgs(varargin,{50,50,...
       prctile(var1(ind),[1,99]),prctile(var2(ind),[1,99]),'','',...
       Trial.stc.list_state_attrib('key'),'default'});


if ischar(states), states = mat2cell(states,1,ones(size(states))); end


%% Plot JPDF(v1,v2), distribution 
hfig = figure(2323);

ind = Trial.stc{'a'};
states(cellfun(@strcmp,states,repmat({'a'},size(states))))=[];
nsts = numel(states);

v1edgs = linspace(v1lim(1),v1lim(2),v1nbins);
v2edgs = linspace(v2lim(1),v2lim(2),v2nbins);

caxmax = [];
N = hist2([var1(ind),var2(ind)],v1edgs,v2edgs);
axes('outerposition',[0,0,1/(nsts+1),1]);
imagesc(v1edgs,v2edgs,N',[0,mean(N(nniz(N(:))))*5]),
axis xy;
ylabel(v2Label);
xlabel(v1Label);
title(ind.label);
xlim(v1edgs([1,end]))
ylim(v2edgs([1,end]))
caxmax(1) = subsref(caxis,substruct('()',{[1],[2]}));

slab = {};
for i = 1:nsts
    hax = axes('outerposition',[i/(nsts+1),0,1/(nsts+1),1]);
    hold on;
    ind = Trial.stc{states{i}};
    slab{end+1} = ind.label; 
    N = hist2([var1(ind),var2(ind)],v1edgs,v2edgs);
    imagesc(v1edgs,v2edgs,N',[0,mean(N(nniz(N(:))))*5]);
    axis xy;
    ylabel(v2Label);
    xlabel(v1Label);
    title(ind.label);
    xlim(v1edgs([1,end]))
    ylim(v2edgs([1,end]))
    caxmax(i+1) = subsref(caxis,substruct('()',{[1],[2]}));
end


ForAllSubplots(['caxis([0,' num2str(median(caxmax)) '])']);

imtype = 'png';

imgfilebase = fullfile(Trial.path.data,'figures',mfilename);
if ~exist(imgfilebase,'dir'), mkdir(imgfilebase),end

%keyboard

% $$$ set(hfig,'Units','inches')
% $$$ set(hfig,'position',get(hfig,'position').*[1,1,0,1]+[0,0,3,0]);
% $$$ set(hfig,'Units','normalized')

%keyboard

set(hax,'Units','inches')
wh = subsref(get(hax,'position'),substruct('()',{[3,4]}));
set(hax,'Units','normalized')
set(hfig,'position',get(hfig,'position').*[0,0,wh(2)/wh(1),1])
pause(.1)

set(hax,'Units','inches')
wh = subsref(get(hax,'position'),substruct('()',{[3,4]}));
set(hax,'Units','normalized')
set(hfig,'position',get(hfig,'position').*[0,0,wh(2)/wh(1),1])


set(hfig,'Units','inches')
set(hfig,'paperposition',get(hfig,'position').*[0,0,1,1]/2)


saveas(hfig,fullfile(Trial.path.data,'figures',mfilename,...
                     [Trial.filebase '.' mfilename '.' Trial.stc.mode '.' tag '.' imtype]),imtype);

imtype = 'eps';
saveas(hfig,fullfile(Trial.path.data,'figures',mfilename,...
                     [Trial.filebase '.' mfilename '.' Trial.stc.mode '.' tag '.' imtype]),'epsc2');


% $$$ savefig(fullfile(Trial.path.data,'figures',mfilename,...
% $$$                  [Trial.filebase '.' mfilename '.' tag '.' imtype]),...
% $$$         hfig,imtype);


%reportfig(fullfile(Trial.path.data,'figures'),figh,'jpdf_ang_height',[],Trial.filebase,100)


