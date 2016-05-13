function req20160310_8_genOptfigs(Trial)
Trial = MTATrial.validate(Trial);


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

% $$$ 
% $$$ figure
% $$$ for s = 1:5;
% $$$ subplot2(3,5,1,s);
% $$$ hold on,
% $$$ plot(acc_ori(:,s).*100)
% $$$ plot(acc_opt(:,s).*100)
% $$$ 
% $$$ subplot2(3,5,2,s);
% $$$ hold on,
% $$$ plot(sen_ori(:,s))
% $$$ plot(sen_opt(:,s))
% $$$ 
% $$$ subplot2(3,5,3,s);
% $$$ hold on,
% $$$ plot(pre_ori(:,s))
% $$$ plot(pre_opt(:,s))
% $$$ end
% $$$ ForAllSubplots('ylim([75,100])')

%rear  opt 1:5
%sit   ori 1:5
%groom opt 1:6
%pause opt 1:5
%walk  opt 1:3

bfets = {bs.bFetInds{1}(1:5);sbind(1:5,2);bs.bFetInds{3}(1:6);bs.bFetInds{4}(1:5);bs.bFetInds{5}(1:3)};

save(fullfile(Trial.spath,[mfilename,'.mat']),'bfets','-v7.3');

%ufets = unique(cell2mat(bfets));



% $$$ for s = 1:numel(bfets),
% $$$     [U,S,V] = svd(dsd.afet(Trial.stc{'a'},bfets{s})'*dsd.afet(Trial.stc{'a'},bfets{s}));
% $$$ figure,
% $$$ out = hist2([dsd.afet(Trial.stc{'a'},bfets{s})*V(:,1),dsd.afet(Trial.stc{'a'},bfets{s})*V(:,2)],edx,edy);
% $$$ edx= 
% $$$ edy= linspace(-2,5,50);
% $$$ out = hist2([dsd.afet(Trial.stc{'a'},bfets{s})*V(:,1),dsd.afet(Trial.stc{'a'},bfets{s})*V(:,2)],edx,edy);
% $$$ caxis([0,40])
% $$$ figure,
% $$$ hist2([dsd.afet(Trial.stc{'r'},bfets{s})*V(:,1),dsd.afet(Trial.stc{'r'},bfets{s})*V(:,2)],edx,edy);
% $$$ caxis([0,40])
% $$$ end

% Rear
% $$$ s = 1;
% $$$ cstate = Trial.stc{'r+s+m+p+w+n&a'};
% $$$ [U,S,V] = svd(dsd.afet(cstate,bfets{s})'*dsd.afet(cstate,bfets{s}));
% $$$ sts = {'rear','s+m+p+w+n&a'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-7,4,60),linspace(-4,3,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(cstate,bfets{s})*V(:,1),dsd.afet(cstate,bfets{s})*V(:,2)],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,1),dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,2)];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[.1,.1],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('PC1')
% $$$ ylabel('PC2')
% $$$ title('Rear Vs All Excluding: Rear')
% $$$ legend({'Rear','All x {Rear}'},'location','northeast');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           false,                               ... Do Not Save FIG
% $$$           'png');                                % Output Format
% $$$ 
% $$$ 
% $$$ 
% $$$ % sit
% $$$ s = 2;
% $$$ cstate = Trial.stc{'s+m+p+w+n&a'};
% $$$ [U,S,V] = svd(dsd.afet(cstate,bfets{s})'*dsd.afet(cstate,bfets{s}));
% $$$ sts = {'s','m+p+w+n&a'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-2.5,2.5,60),linspace(-2.2,3,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(cstate,bfets{s})*V(:,1),dsd.afet(cstate,bfets{s})*V(:,2)],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,1),dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,2)];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('PC1')
% $$$ ylabel('PC2')
% $$$ title('Sit Vs All Excluding: Rear,Sit')
% $$$ legend({'Sit','all x {Rear,Sit}'},'location','northeast');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           false,                               ... Do Not Save FIG
% $$$           'png');                                % Output Format
% $$$ 
% $$$ 
% $$$ 
% $$$ s = 3;
% $$$ cstate = Trial.stc{'m+p+w+n&a'};
% $$$ [U,S,V] = svd(dsd.afet(cstate,bfets{s})'*dsd.afet(cstate,bfets{s}));
% $$$ sts = {'m','p+w+n&a'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-3,3,60),linspace(-2.5,3,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(cstate,bfets{s})*V(:,1),dsd.afet(cstate,bfets{s})*V(:,2)],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,1),dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,2)];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('PC1')
% $$$ ylabel('PC2')
% $$$ title('Groom Vs All excluding: Rear,Sit,Groom')
% $$$ legend({'groom','all x {rear,sit,groom}'},'location','northwest');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           false,                               ... Do Not Save FIG
% $$$           'png');                                % Output Format
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ s = 4;
% $$$ cstate = Trial.stc{'p+w+n&a'};
% $$$ [U,S,V] = svd(dsd.afet(cstate,bfets{s})'*dsd.afet(cstate,bfets{s}));
% $$$ sts = {'p','w+n'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-5.5,3,60),linspace(-2.5,2.5,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(cstate,bfets{s})*V(:,1),dsd.afet(cstate,bfets{s})*V(:,2)],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,1),dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,2)];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('PC1')
% $$$ ylabel('PC2')
% $$$ title('Pause Vs All excluding: Rear,Sit,Groom,Pause')
% $$$ legend({'Pause','All x {Rear,Sit,Groom,Pause}'},'location','northeast');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           false,                               ... Do Not Save FIG
% $$$           'png');                                % Output Format
% $$$ 
% $$$ 
% $$$ 
% $$$ s = 5;
% $$$ [U,S,V] = svd(dsd.afet(Trial.stc{'w+n'},bfets{s})'*dsd.afet(Trial.stc{'w+n'},bfets{s}));
% $$$ sts = {'walk','turn'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-4,1.5,60),linspace(-2.5,2,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(Trial.stc{'w+n'},bfets{s})*V(:,1),dsd.afet(Trial.stc{'w+n'},bfets{s})*V(:,2)],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,1),dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,2)];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('PC1')
% $$$ ylabel('PC2')
% $$$ title('Walk Vs Turn')
% $$$ legend({'Walk','Turn'},'location','northeast');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           false,                               ... Do Not Save FIG
% $$$           'png');                                % Output Format




% $$$ % Rear
% $$$ [U,S,V] = svd(dsd.afet(cstate,bfets{s})'*dsd.afet(cstate,bfets{s}));
% $$$ sts = {'rear','s+m+p+w+n&a'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-7,4,60),linspace(-4,3,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(cstate,bfets{s})*V(:,1),dsd.afet(cstate,bfets{s})*V(:,2)],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,1),dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,2)];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o/sum(o(:))',[.1,.1],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ 
% $$$ 
% $$$ s = 1;
% $$$ cstate = Trial.stc{'r&a'};
% $$$ o = hist2([dsd.afet(cstate,bfets{s})*V(:,1),dsd.afet(cstate,bfets{s})*V(:,2)],hedgs{:});
% $$$ figure,imagesc([o/sum(o(:))]');

% $$$ cstate = Trial.stc{'s+m+p+w+n&a'};
% $$$ o = hist2([dsd.afet(cstate,bfets{s})*V(:,1),dsd.afet(cstate,bfets{s})*V(:,2)],hedgs{:});
% $$$ figure,imagesc([o/sum(o(:))]');
% $$$ 
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,1),dsd.afet(Trial.stc{sts{i}},bfets{s})*V(:,2)];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = o/sum(o(:));
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[.001,.001],'linewidth',2.5,'Color',stc(i))
% $$$ 
