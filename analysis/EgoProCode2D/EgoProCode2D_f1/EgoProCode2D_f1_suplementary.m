



% pairwise difference descending vs trough
figure,
hold('on');
grid('on');
ecdf(tMeanPos(uidsCA1,2) - tMeanPos(uidsCA1,1))
ecdf(tMeanPos(uidsCA1,3) - tMeanPos(uidsCA1,2))
ecdf(tMeanPos(uidsCA1,3) - tMeanPos(uidsCA1,1))
title({'Emperical Cumulative Density Function',...
       'of pairwise difference of anteroposteror',...
       'ordinate between CA1 phase groups'});
legend({'trough - descend', 'ascend - trough', 'ascend - descend'},'Location','southeast');
xlabel('Difference (cm)');
ylabel('Probability');
xlim([-5,20]);

figure,
hold('on');
grid('on');
ecdf(tMeanPos(uidsCA3,2) - tMeanPos(uidsCA3,1));
ecdf(tMeanPos(uidsCA3,3) - tMeanPos(uidsCA3,2));
ecdf(tMeanPos(uidsCA3,3) - tMeanPos(uidsCA3,1));
title({'Emperical Cumulative Density Function',...
       'of pairwise difference of anteroposteror',...
       'ordinate between CA3 phase groups'});
legend({'trough - descend', 'ascend - trough', 'ascend - descend'},'Location','southeast');
xlabel('Difference (cm)');
ylabel('Probability');
xlim([-5,20]);