function [occMap,Bins] = xyocc(Trial,varargin)
% function occ = xyocc(Trial,pos,type,binDims)
% [pos,type,binDims] = DefaultArgs(varargin,{Trial.stc{'w'},'xy',[20,20]
% Always us xy for now
[tpos,type,binDims,SmoothingWeights,bound_lims,display,report] = DefaultArgs(varargin,{Trial.stc.states,'xy',[20,20],[1.8,1.8],[],true,true});


if isempty(tpos),
    Trial = labelAllBhv(Trial);
    tpos = Trial.stc.states;
end


if ~iscell(tpos)
    tpos = {tpos};
end

for t = 1:numel(tpos),
    
    if isa(tpos{t},'MTADepoch');       
        xyz = Trial.load('xyz');
        posSampleRate = xyz.sampleRate;
        pos = sq(xyz(tpos{t},Trial.trackingMarker,ismember('xyz',type)));
    else
        pos = tpos{t};
    end

    ndims = length(type);

    if isempty(bound_lims),
        bound_lims = Trial.maze.boundaries(ismember('xyz',type),:);
    end

    Nbin = round(abs(diff(bound_lims,1,2))./binDims');


    Bins = cell(1,ndims);
    k = Nbin./abs(diff(bound_lims,1,2));
    msize = round(abs(diff(bound_lims,1,2)).*k);
    for i = 1:ndims
        Bins{i} = ([1:msize(i)]'-1)./repmat(k(i)',msize(i),1)+repmat(bound_lims(i,1),msize(i),1)+round(repmat(k(i)',msize(i),1).^-1/2);
    end


    Pos = round((pos-repmat(bound_lims(:,1)',size(pos,1),1)).*repmat(k',size(pos,1),1))+1;
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

    OccThresh = 0.06;%0.12;%
    gtind = occ>OccThresh;

    occm = NaN(msize');
    occm(gtind) = occ(gtind);
    occm(~gtind) = NaN;

    occMap{t} = occm;

end

if display
    
    hfig = figure(48482),clf
    y = ceil(numel(tpos)/4);
    if numel(tpos)>=4, 
        x = 4;
    else
        x = numel(tpos);
    end
    
    for t = 1:numel(tpos),
        subplot(y,x,t);

        imagescnan(cat(2,Bins,occMap{t}'),[],[],true,[0,0,0]);
        if isa(tpos{t},'MTADepoch'),
        title([Trial.filebase,':',tpos{t}.label])
        else
            title([Trial.filebase])
        end
    end
 
    set(hfig,'position',get(hfig,'position').*[1,1,0,0]+[0,0,x*330,y*300])
    
    if report
        reportfig(fullfile(Trial.path.data,'figures'),hfig,['xyocc'],false,Trial.filebase,100);
    end

end

    