function [occMap,Bins] = generate_occupancy_map(Trial,varargin)
% function occ = xyocc(Trial,pos,type,binDims)
% [pos,type,binDims] = DefaultArgs(varargin,{Trial.stc{'w'},'xy',[20,20]
% Always us xy for now

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('pos',              Trial.load('xyz'),                                          ...
                 'state',            [],                                                         ...
                 'type',             'xy',                                                       ...
                 'binDims',          [20,20],                                                    ...
                 'SmoothingWeights', [2.2,2.2],                                                  ...
                 'boundLims',        [],                                                         ...
                 'OccThresh',        [0.06],                                                     ...
                 'display',          false,                                                      ...
                 'report',           false                                                       ...
);
[pos,state,type,binDims,SmoothingWeights,boundLims,OccThresh,display,report] = ...
    DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
if isempty(state),
    Trial = labelAllBhv(Trial);
    state = Trial.stc.states;
end


if ~iscell(state)
    state = {state};
end

for t = 1:numel(state),
    
    posSampleRate = pos.sampleRate;
    pos = sq(pos(state{t},Trial.trackingMarker,ismember('xyz',type)));

    ndims = length(type);

    if isempty(boundLims),
        boundLims = Trial.maze.boundaries(ismember('xyz',type),:);
    end

    Nbin = round(abs(diff(boundLims,1,2))./binDims');


    Bins = cell(1,ndims);
    k = Nbin./abs(diff(boundLims,1,2));
    msize = round(abs(diff(boundLims,1,2)).*k);
    for i = 1:ndims,
        Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+...
                  repmat(boundLims(i,1),msize(i),1)+...
                  round(repmat(k(i)',msize(i),1).^-1/2);
    end


    Pos = round((pos-repmat(boundLims(:,1)',size(pos,1),1)).*repmat(k',size(pos,1),1))+1;
    for i = 1:ndims   
        Pos(Pos(:,i)<1|Pos(:,i)>Nbin(i)|~nniz(Pos),:) = [];
    end
    occ = accumarray(Pos,1,msize')./posSampleRate;


    sind = cell(1,ndims);
    for i = 1:ndims,
        sind{i} = linspace(-round(msize(i)/2),round(msize(i)/2),msize(i));
    end
    [sind{:}] = ndgrid(sind{:});
    for i = 1:ndims,
        sind{i} = sind{i}.^2/SmoothingWeights(i)^2/2;
    end
    Smoother = exp(sum(-cat(ndims+1,sind{:}),ndims+1));
    Smoother = Smoother./sum(Smoother(:));
    occ = convn(occ,Smoother,'same');

    %OccThresh = 0.06;%0.12;%
    gtind = occ>OccThresh;

    occm = NaN(msize');
    occm(gtind) = occ(gtind);
    occm(~gtind) = NaN;

    occMap{t} = occm;

end

if display
    
    hfig = figure(48482),clf
    y = ceil(numel(state)/4);
    if numel(state)>=4, 
        x = 4;
    else
        x = numel(state);
    end
    
    for t = 1:numel(state),
        subplot(y,x,t);

        imagescnan(cat(2,Bins,occMap{t}'),[],[],true,[0,0,0]);
        if isa(state{t},'MTADepoch'),
        title([Trial.filebase,':',state{t}.label])
        else
            title([Trial.filebase])
        end
    end
 
    set(hfig,'position',get(hfig,'position').*[1,1,0,0]+[0,0,x*330,y*300])
    
    if report
        reportfig(fullfile(Trial.path.data,'figures'),hfig,['xyocc'],false,Trial.filebase,100);
    end

end

% END MAIN -----------------------------------------------------------------------------------------
