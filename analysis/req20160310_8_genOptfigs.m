function req20160310_8_genOptfigs(Trial)
Trial = MTATrial.validate(Trial);
mfilename = 'req20160310_8_genOptfigs';

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


% $$$ figure
% $$$ set(hfig,'PaperPositionMode','auto');
% $$$ set(gcf,'units','centimeters');
% $$$ set(gcf,'Position',[0,0,25,5.5])
% $$$ set(gcf,'PaperUnits','centimeters');
% $$$ set(gcf,'PaperPosition',[0,0,25,5.5])
% $$$ 
% $$$ for s = 1:5;
% $$$ subplot2(1,5,1,s);
% $$$ hold on,
% $$$ %acc
% $$$ hs = plot(acc_ori(:,s).*100);
% $$$ hs.Color = [1,0.75,0.75];
% $$$ hs = plot(acc_opt(:,s).*100);
% $$$ hs.Color = [1,0,0];
% $$$ 
% $$$ % $$$ %sen
% $$$ % $$$ hs = plot(sen_ori(:,s));
% $$$ % $$$ hs.Color = [0.75,0.75,1];
% $$$ % $$$ hs = plot(sen_opt(:,s));
% $$$ % $$$ hs.Color = [0,0,1];
% $$$ % $$$ 
% $$$ % $$$ %pre
% $$$ % $$$ hs = plot(pre_ori(:,s));
% $$$ % $$$ hs.Color = [0.75,1,0.75];
% $$$ % $$$ hs = plot(pre_opt(:,s));
% $$$ % $$$ hs.Color = [0,1,0];
% $$$ title(dsd.stateOrd{s})
% $$$ yl = ylim;
% $$$ ylim([yl(1),100]);
% $$$ xlim([1,30])
% $$$ set(gca,'units','centimeters');
% $$$ cpos = get(gca,'Position');
% $$$ set(gca,'Position',[cpos(1),cpos(2),2.5,3.5])
% $$$ end
% $$$     set(findall(gcf,'-property','FontSize'),'FontSize',8)
% $$$ print(gcf,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_3',...
% $$$                      'fig3_B_accuracy.eps'))


% jg05-20120317
%rear  opt 1:5
%sit   ori 1:5
%groom opt 1:6
%pause opt 1:5
%walk  opt 1:3
bfets = {bs.bFetInds{1}(1:5);sbind(1:5,2);bs.bFetInds{3}(1:6);bs.bFetInds{4}(1:5);bs.bFetInds{5}(1:3)};

% Ed03-20140625
%bfets = {bs.bFetInds{1}(1:7);bs.bFetInds{1}(1:2);bs.bFetInds{3}(1:17);bs.bFetInds{4}(1:14);bs.bFetInds{5}(1:9)};


save(fullfile(Trial.spath,[mfilename,'.mat']),'bfets','-v7.3');

%ufets = unique(cell2mat(bfets));


% $$$ fets_jg = [2 4 6 7 8 39 40 41 42 43 45 46 47 53 56 58 59]
% $$$ 
% $$$ fets_Ed = [1 2 3 4 5 6 8 10 27 30 35 37 38 39 40 41 42 43 44 45 46 ...
% $$$            47 48 49 50 57 58 59];



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



% $$$ % Rear
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
% $$$           true,                                ... Save FIG
% $$$           'png');                                % Output Format
% $$$ 

% $$$ s = 1;
% $$$ cstate = Trial.stc{'r+s+m+p+w+n&a'};
% $$$ sts = {'rear','s+m+p+w+n&a'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-2,4.5,60),linspace(-2.5,2.5,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ set(hfig,'units','centimeters');
% $$$ set(hfig,'Position',[0,0,8,8])
% $$$ hist2([dsd.afet(cstate,bfets{s}(1)),dsd.afet(cstate,bfets{s}(3))],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(1)),dsd.afet(Trial.stc{sts{i}},bfets{s}(3))];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[5,5],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('Normalized Spine Pitch (A.U.)')
% $$$ ylabel('Normalized Head Speed Z-axis (A.U.)')
% $$$ title('Rear Vs All Excluding: Rear')
% $$$ % $$$ legend({'Rear','All x {Rear}'},'location','northeast');
% $$$ % $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$ % $$$           hfig,                                ... Figure handle
% $$$ % $$$           [mfilename,'-fet'],                         ... Figure Set Name
% $$$ % $$$           'req',                               ... Directory where figures reside
% $$$ % $$$           false,                               ... Do Not Preview
% $$$ % $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$ % $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$ % $$$           [],                                  ... Resolution
% $$$ % $$$           true,                                ... Save FIG
% $$$ % $$$           'png');                                % Output Format
% $$$ set(gca,'units','centimeters');
% $$$ set(gca,'Position',[4-3.5/2,4-3.5/2,3.5,3.5])
% $$$ set(findall(gcf,'-property','FontSize'),'FontSize',8)
% $$$ print(hfig,'-depsc2',fullfile(getenv('PROJECT'),'manuscripts/man2015-jgEd-MoCap/Figures/Figure_2',...
% $$$                      'fig3_C_jpdf_rear.eps'))
          

% $$$ 
% $$$ % sit
% $$$ % PCA
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
% $$$ % Actual features
% $$$ s = 2;
% $$$ cstate = Trial.stc{'s+m+p+w+n&a'};
% $$$ sts = {'s','m+p+w+n&a'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-2.5,3,60),linspace(-2.2,2,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(cstate,bfets{s}(1)),dsd.afet(cstate,bfets{s}(2))],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(1)),dsd.afet(Trial.stc{sts{i}},bfets{s}(2))];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[5,5],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('Normalized Height of Lower Body (A.U.)')
% $$$ ylabel('Normalized Height of Pelvis (A.U.)')
% $$$ title('Sit Vs All Excluding: Rear,Sit')
% $$$ legend({'Sit','all x {Rear,Sit}'},'location','northeast');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename,'-fet'],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           true,                                ... Save FIG
% $$$           'png');                                % Output Format
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
% $$$ s = 3;
% $$$ cstate = Trial.stc{'m+p+w+n&a'};
% $$$ sts = {'m','p+w+n&a'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-3,4,60),linspace(-1.5,4,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(cstate,bfets{s}(1)),dsd.afet(cstate,bfets{s}(2))],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(1)),dsd.afet(Trial.stc{sts{i}},bfets{s}(2))];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('Normalized Height of Lower Body (A.U.)')
% $$$ ylabel('Normalized Spine Curvature in XY Plane (A.U.)')
% $$$ title('Groom Vs All excluding: Rear,Sit,Groom')
% $$$ legend({'groom','all x {rear,sit,groom}'},'location','northwest');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename,'-fet'],                         ... Figure Set Name
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
% $$$ s = 4;
% $$$ cstate = Trial.stc{'p+w+n&a'};
% $$$ sts = {'p','w+n'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-2,2.5,60),linspace(-2,3.5,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(cstate,bfets{s}(2)),dsd.afet(cstate,bfets{s}(5))],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(2)),dsd.afet(Trial.stc{sts{i}},bfets{s}(5))];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('Normalized Lower Body Height (A.U.)')
% $$$ ylabel('Normalized Lower Body Speed (A.U.)')
% $$$ title('Pause Vs All excluding: Rear,Sit,Groom,Pause')
% $$$ legend({'Pause','All x {Rear,Sit,Groom,Pause}'},'location','northeast');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename,'-fet'],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           true,                                ... Do Not Save FIG
% $$$           'png');                                % Output Format
          
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
% $$$ 
% $$$ s = 5;
% $$$ sts = {'walk','turn'};
% $$$ stc = 'rc';
% $$$ hedgs = {linspace(-1,2.5,60),linspace(-1.5,3.5,60)};
% $$$ edgs  = hedgs;
% $$$ [edgs{:}] = get_histBinCenters(edgs);
% $$$ [X,Y] = meshgrid(edgs{:});
% $$$ hfig = figure(201603108);clf;
% $$$ hist2([dsd.afet(Trial.stc{'w+n'},bfets{s}(2)),dsd.afet(Trial.stc{'w+n'},bfets{s}(3))],hedgs{:});
% $$$ caxis([0,40])
% $$$ hold on,
% $$$ for i = 1:numel(sts),
% $$$     b = [dsd.afet(Trial.stc{sts{i}},bfets{s}(2)),dsd.afet(Trial.stc{sts{i}},bfets{s}(3))];
% $$$     o = hist2(b,hedgs{:});
% $$$     F = [.05 .1 .05; .1 .4 .1; .05 .1 .05];
% $$$     o = conv2(o,F,'same');
% $$$     contour(X,Y,o',[3,3],'linewidth',2.5,'Color',stc(i))
% $$$ end
% $$$ xlabel('Normalized Head Speed (A.U.)')
% $$$ ylabel('Normalized All Trajectories PPC (A.U.)')
% $$$ title('Walk Vs Turn')
% $$$ legend({'Walk','Turn'},'location','northeast');
% $$$ 
% $$$ reportfig(fullfile(getenv('PROJECT'),'figures'),  ... Path where figures are stored
% $$$           hfig,                                ... Figure handle
% $$$           [mfilename,'-fet'],                         ... Figure Set Name
% $$$           'req',                               ... Directory where figures reside
% $$$           false,                               ... Do Not Preview
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Tumbnail caption
% $$$           [Trial.filebase '-' dsd.stateOrd{s}],... Expanded caption
% $$$           [],                                  ... Resolution
% $$$           true,                                ... Do Not Save FIG
% $$$           'png');                                % Output Format
% $$$           



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
