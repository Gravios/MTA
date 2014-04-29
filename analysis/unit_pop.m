sesList = {{'jg05-20120309','cof','all'},...
           {'jg05-20120310','cof','all'},...
           {'jg05-20120317','cof','all'}};



for ses = 1:numel(sesList),
    Trial = MTATrial(sesList{ses}{1},sesList{ses}{3},sesList{ses}{2});
    Trial.xyz.load(Trial);
    Trial.load('nq');
    nq{ses} = Trial.nq;
end

units = find(Trial.nq.SpkWidthR>0.8&Trial.nq.eDist>18)';

figure,plot(anq.SpkWidthR(anq.eDist>30),anq.bRat(anq.eDist>30),'.')
figure,plot(anq.SpkWidthR(anq.eDist>30),anq.AmpSym(anq.eDist>30),'.')

