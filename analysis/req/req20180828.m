
% Stripping the kloosterman bayesian decoding 
% stats on two line replacement 

dbstop at 176 in bhv_decode.m

t = zeros([1000,1]);
bins = ones([1,bufferSize]);
tufr = ufr(:,:)-eps;
%rateMap = rateMap+1e-3;
for n = 1:1000;
    tic;
    E = exp(-sum(rateMap,2)*bins*spikeWindow+log(rateMap)*(ufr.data(ind,:))');
    E = bsxfun(@rdivide,E,sum(E));
    t(n) = toc;
end

fE = mazeMask;
fE = nan(size(mazeMask));
fE(mazeMask(:)) = E(:,400);
fE(mazeMask(:)) = tE(:,400);
figure,
imagesc(fE(:,:)');

figure;plot(E(:,200),'LineWidth',5);hold('on');plot(tE(:,200),'LineWidth',2);

o = zeros([1000,1]);
for n = 1:1000;
    tic
    tE = decode_bayesian_poisson(rateMap,...
                                 (ufr.data(ind,:)'-eps),...
                                 'alpha',1,...
                                 'bins',ones([1,bufferSize])*spikeWindow,...
                                 'baseline',1e-3);
    o(n) = toc;
end

figure();
hold('on');
eds = linspace(0,0.3,100);
hax = bar(eds,histc(o,eds),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.4;
hax.EdgeAlpha = 0.4;
hax = bar(eds,histc(t,eds),'histc');
hax.FaceColor = 'c';
hax.EdgeColor = 'c';
hax.FaceAlpha = 0.4;
hax.EdgeAlpha = 0.4;
