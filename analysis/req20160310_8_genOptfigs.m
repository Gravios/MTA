function req20160310_8_genOptfigs(Trial)
Trial = MTATrial.validate


acc_ori = [];
acc_opt = [];

sen_ori = [];
sen_opt = [];

pre_ori = [];
pre_opt = [];

sbind = [];

dsd = load(fullfile(Trial.spath,'req20160310_1_preproc-afet.mat'));
bs = load(fullfile(Trial.spath,'req20160310_5_genfigs.mat'));
for s = 1:5
       
    oind = [repmat([1:59],1,2)',zeros([118,1])];
    aind = oind(:,1);
    for sh = 1:117,
        oind = [oind;[circshift(aind,-sh),aind]];
    end
    slind = oind(dsd.fetInds{s},:);
    ofet =reshape(slind,[],1);
    best_inds = histc(ofet,1:59);
    [~,sbind(:,s)] = sort(best_inds,'descend');
    
    ori = load(fullfile(Trial.spath,['req20160310_4_accumStats',num2str(s),'.mat']));    
    opt = load(fullfile(Trial.spath,['req20160310_7_accumOptStats',num2str(s),'.mat']));

    acc_ori(:,s) =  cell2mat({ori.accum_acc});
    acc_opt(:,s) =  cell2mat({opt.accum_acc});    

    sen_ori(:,s) =  cell2mat({ori.accum_sen});
    sen_opt(:,s) =  cell2mat({opt.accum_sen});    

    pre_ori(:,s) =  cell2mat({ori.accum_pre});
    pre_opt(:,s) =  cell2mat({opt.accum_pre});    

end


figure
for s = 1:5;
subplot2(3,5,1,s);
hold on,
plot(acc_ori(:,s).*100)
plot(acc_opt(:,s).*100)

subplot2(3,5,2,s);
hold on,
plot(sen_ori(:,s))
plot(sen_opt(:,s))

subplot2(3,5,3,s);
hold on,
plot(pre_ori(:,s))
plot(pre_opt(:,s))
end
ForAllSubplots('ylim([75,100])')

%rear  opt 1:5
%sit   ori 1:5
%groom opt 1:6
%pause opt 1:5
%walk  opt 1:3
bfets = [bs.bFetInds{1}(1:5);sbind(1:5,2);bs.bFetInds{3}(1:6);bs.bFetInds{4}(1:5);bs.bFetInds{5}(1:9)];
bfets = unique(bfets);


