function [xcomp,ycomp,zcomp] = compute_2d_discrete_func(xcomp,ycomp,vals,funcHandle)

xcomp.edgs = xcomp.edgs(:)';
ycomp.edgs = ycomp.edgs(:)';

xcomp.ctrs = mean([xcomp.edgs(1:end-1);xcomp.edgs(2:end)]);
ycomp.ctrs = mean([ycomp.edgs(1:end-1);ycomp.edgs(2:end)]);

if ~isfield(xcomp,'inds')
    xcomp.inds = discretize(xcomp.data,xcomp.edgs);
end
if ~isfield(ycomp,'inds')
    ycomp.inds = discretize(ycomp.data,ycomp.edgs);
end

zcomp.xctr = xcomp.ctrs;
zcomp.yctr = ycomp.ctrs;
nind = nniz(xcomp.inds) & nniz(ycomp.inds) & nniz(vals);
zcomp.data = accumarray([xcomp.inds(nind),ycomp.inds(nind)],vals(nind),[numel(xcomp.ctrs),numel(xcomp.ctrs)],funcHandle);

zcomp.count = accumarray([xcomp.inds(nind),ycomp.inds(nind)],ones([sum(nind),1]),[numel(xcomp.ctrs),numel(xcomp.ctrs)],@sum);