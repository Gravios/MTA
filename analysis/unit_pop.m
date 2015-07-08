


sesList = {{'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};



for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
    Trial.xyz.load(Trial);
    Trial.load('nq');
    nq{ses} = Trial.nq;
end
nq = cat(1,nq{:});
anq = CatStruct(nq);
units = find(Trial.nq.SpkWidthR>0.8&Trial.nq.eDist>18)';



figH = figure(48849),
plot(anq.SpkWidthR(anq.eDist>15),anq.AmpSym(anq.eDist>15),'.')
title('Fit for hyperplane perpendicular to selected dimensions)
xlabel('SpkWidthR');
ylabel('AmpSym');
pram = draw_lines(figH,'line_fit');
saveas(figH,'/gpfs01/sirota/bach/homes/gravio/figures/Unit_Selection/spkWR_AmpSym.png');
saveas(figH,'/gpfs01/sirota/bach/homes/gravio/figures/Unit_Selection/spkWR_AmpSym.fig');
