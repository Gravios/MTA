m

figure,
for t = 1:30
    subplot(5,6,t);
    %vxy = fet_href_HXY(Trials{t});
    vxy = fet_HB_pitch(Trials{t},30);
    out{t} = hist2([vxy([Trials{t}.stc{'r'}],3),vxy([Trials{t}.stc{'r'}],1)],...
                   linspace([-1.5,1.5,24]),...
                   linspace([-0.5,2,24]));
end

% $$$     hist2([vxy([Trials{t}.stc{'r'}],2),vxy([Trials{t}.stc{'r'}],1)],...
% $$$           linspace([-50,50,24]),...
% $$$           linspace([-25,75,24]));

bin.headPitch.edges = linspace([-1.5,1.5,24]);
bin.headPitch.centers = mean([bin.headPitch.edges(1:end-1);bin.headPitch.edges(2:end)]);
bin.headPitch.count = numel(bin.headPitch.centers);
bin.headPitch.label = 'Head Pitch HC-NS (rad)';

bin.headHeight.edges = linspace([0,300,24]);
bin.headHeight.centers = mean([bin.headHeight.edges(1:end-1);bin.headHeight.edges(2:end)]);
bin.headHeight.count = numel(bin.headHeight.centers);
bin.headHeight.label = 'Head Height HC (mm)';

bin.bodyPitch.edges = linspace([-0.5,2,24]);
bin.bodyPitch.centers = mean([bin.bodyPitch.edges(1:end-1);bin.bodyPitch.edges(2:end)]);
bin.bodyPitch.count = numel(bin.bodyPitch.centers);
bin.bodyPitch.label = 'Body Pitch BP-BU (rad)';


dlabels = {'interior','perimeter'};
dsymbols ={'<','>'};
dthresh = 300;

% $$$ fet = cf(@(T) fet_HB_pitch(T,30), Trials);
% $$$ mxyz = cf(@(x) x.copy(),xyz);
% $$$ cf(@(x) x.resample(30), mxyz);

% headPitch vs headHeight
jpdf.headPitch_headHeight.out = {};
dlogic = {@lt,@gt}
for t = 1:30
    for d = 1:2
    disp(['[INFO] processing ',Trials{t}.filebase,' ...'])
    dist = sqrt(sum(mxyz{t}(:,'hcom',[1,2]).^2,3));
    rper = cast(Trials{t}.stc{'r'},'TimeSeries');
    rper.resample(mxyz{t});
    ind = rper.data & dlogic{d}(dist,dthresh);
    jpdf.headPitch_headHeight.out{t}(:,:,d) =  ...
        hist2([fet{t}(ind,3),mxyz{t}(ind,'hcom',3)],...
               bin.headPitch.edges,                                      ...
               bin.headHeight.edges                                       ...
             );
    end
end


for d = 1:2;
jname = 'headPitch_headHeight'; 
jpdf.(jname).x.fname = 'headPitch';   jpdf.(jname).x.lname = 'head pitch';
jpdf.(jname).y.fname = 'headHeight';  jpdf.(jname).y.lname = 'head height';
jpdf.(jname).bin.(jpdf.(jname).x.fname) = bin.(jpdf.(jname).x.fname);
jpdf.(jname).bin.(jpdf.(jname).y.fname) = bin.(jpdf.(jname).y.fname);
jpdf.(jname).data = sum(cat(4,jpdf.(jname).out{:}),4);
jpdf.(jname).label = {'Joint Probability Distribution', ...
                    ['between ',jpdf.(jname).x.lname,' and ',jpdf.(jname).y.lname],...
                    'during rearing in the arena perimeter',...
                    ['head distance from center ',dsymbols{d},' ',num2str(dthresh),'mm, maze radius: 400mm']};
figure();
imagesc(jpdf.(jname).bin.(jpdf.(jname).x.fname).centers,    ...
        jpdf.(jname).bin.(jpdf.(jname).y.fname).centers,    ...
        jpdf.(jname).data(:,:,d)'                                  ...
        );
xlabel(bin.(jpdf.(jname).x.fname).label);
ylabel(bin.(jpdf.(jname).y.fname).label);
title(jpdf.(jname).label);
axis(gca(),'xy');
colormap(gca(),'jet');
set(gca(),'ColorScale','log');
cax = colorbar();
ylabel(cax,'Count');
saveas(gcf(),fullfile('/storage/gravio/figures/analysis/afig/',['jpdf_',jname,'_REAR_',dlabels{d},'.pdf']));
end



out = {};
iut = {};
for t = 3:30
    disp(['[INFO] processing ',Trials{t}.filebase,' ...'])
    dist = sqrt(sum(mxyz{t}(:,'hcom',[1,2]).^2,3));
    rper = cast(Trials{t}.stc{'r'},'TimeSeries');
    rper.resample(mxyz{t});
    ind = rper.data & dist<300;
    iut{t} =  hist2([mxyz{t}(ind,'hcom',3),fet{t}(ind,1)],...
                   bin.headHeight.edges,                                      ...
                   bin.bodyPitch.edges                                       ...
                   );
    ind = rper.data & dist>300;
    out{t} =  hist2([mxyz{t}(ind,'hcom',3),fet{t}(ind,1)],...
                   bin.headHeight.edges,                                      ...
                   bin.bodyPitch.edges                                       ...
                   );
end


jpdf.headPitch_headHeight.bin.headPitch = bin.headPitch;
jpdf.headPitch_headHeight.bin.headHeight = bin.headHeight;
jpdf.headPitch_headHeight.data = sum(cat(3,out{:}),3);
jpdf.headPitch_headHeight.label = {'Joint Probability Distribution', ...
                    'between head pitch and head height',...
                   'during rearing'};

figure();

imagesc(jpdf.headPitch_headHeight.bin.headPitch.centers,    ...
        jpdf.headPitch_headHeight.bin.headHeight.centers,    ...
        jpdf.headPitch_headHeight.data'                       );

xlabel(bin.headPitch.label);
ylabel(bin.headHeight.label);
title(jpdf.headPitch_headHight.label);

axis(gca(),'xy');
colormap(gca(),'jet');
set(gca(),'ColorScale','log');
cax = colorbar();
ylabel(cax,'Count');





jpdf.headPitch_headHeight.bin.headPitch = bin.headPitch;
jpdf.headPitch_headHeight.bin.headHeight = bin.headHeight;
jpdf.headPitch_headHeight.data = sum(cat(3,iut{:}),3);
jpdf.headPitch_headHeight.label = {'Joint Probability Distribution', ...
                    'between head pitch and head height',...
                    'during rearing in the arena perimeter',...
                    'head distance from center <30cm, maze radius: 40cm'};
figure();
imagesc(jpdf.headPitch_headHeight.bin.headPitch.centers,    ...
        jpdf.headPitch_headHeight.bin.headHeight.centers,    ...
        jpdf.headPitch_headHeight.data'                       );
xlabel(bin.headPitch.label);
ylabel(bin.headHeight.label);
title(jpdf.headPitch_headHeight.label);
axis(gca(),'xy');
colormap(gca(),'jet');
set(gca(),'ColorScale','log');
cax = colorbar();
ylabel(cax,'Count');




jpdf.headPitch_headHeight.bin.headPitch = bin.headPitch;
jpdf.headPitch_headHeight.bin.headHeight = bin.headHeight;
jpdf.headPitch_headHeight.data = sum(cat(3,out{:}),3);
jpdf.headPitch_headHeight.label = {'Joint Probability Distribution', ...
                    'between head pitch and head height',...
                    'during rearing in the arena perimeter',...
                    'head distance from center >30cm, maze radius: 40cm'};

figure();
imagesc(jpdf.headPitch_headHeight.bin.headPitch.centers,    ...
        jpdf.headPitch_headHeight.bin.headHeight.centers,    ...
        jpdf.headPitch_headHeight.data'                       );
xlabel(bin.headPitch.label);
ylabel(bin.headHeight.label);
title(jpdf.headPitch_headHeight.label);
axis(gca(),'xy');
colormap(gca(),'jet');
set(gca(),'ColorScale','log');
cax = colorbar();
ylabel(cax,'Count');





jpdf.headPitch_bodyPitch.bin.headPitch = bin.headPitch;
jpdf.headPitch_bodyPitch.bin.bodyPitch = bin.bodyPitch;
jpdf.headPitch_bodyPitch.data = sum(cat(3,out{:}),3);
jpdf.headPitch_bodyPitch.label = {'Joint Probability Distribution', ...
                    'between head pitch and body pitch',...
                   'during rearing'};

figure();

imagesc(jpdf.headPitch_bodyPitch.bin.headPitch.centers,    ...
        jpdf.headPitch_bodyPitch.bin.bodyPitch.centers,    ...
        jpdf.headPitch_bodyPitch.data'                       );

xlabel(bin.headPitch.label);
ylabel(bin.bodyPitch.label);
title(jpdf.headPitch_bodyPitch.label);

axis(gca(),'xy');
colormap(gca(),'jet');
set(gca(),'ColorScale','log');
cax = colorbar();
ylabel(cax,'Count');



jname = 'headHeight_bodyPitch'; d = 2;
jpdf.(jname).x.fname = 'headHeight'; jpdf.(jname).x.lname = 'head height';
jpdf.(jname).y.fname = 'bodyPitch';  jpdf.(jname).y.lname = 'body pitch';
jpdf.(jname).bin.(jpdf.(jname).x.fname) = bin.(jpdf.(jname).x.fname);
jpdf.(jname).bin.(jpdf.(jname).y.fname) = bin.(jpdf.(jname).y.fname);
jpdf.(jname).data = sum(cat(4,out{:}),4);
jpdf.(jname).label = {'Joint Probability Distribution', ...
                    ['between ',jpdf.(jname).x.lname,' and ',jpdf.(jname).y.lname],...
                    'during rearing in the arena perimeter',...
                    'head distance from center >30cm, maze radius: 40cm'};
figure();
imagesc(jpdf.(jname).bin.(jpdf.(jname).x.fname).centers,    ...
        jpdf.(jname).bin.(jpdf.(jname).y.fname).centers,    ...
        jpdf.(jname).data(:,:,d)'                                  ...
        );
xlabel(bin.(jpdf.(jname).x.fname).label);
ylabel(bin.(jpdf.(jname).y.fname).label);
title(jpdf.(jname).label);
axis(gca(),'xy');
colormap(gca(),'jet');
set(gca(),'ColorScale','log');
cax = colorbar();
ylabel(cax,'Count');
saveas(gcf(),fullfile('/storage/gravio/figures/analysis/afig/',['jpdf_',jname,'_REAR_perimeter.pdf']));



jname = 'headHeight_bodyPitch';
jpdf.(jname).x.fname = 'headHeight'; jpdf.(jname).x.lname = 'head height';
jpdf.(jname).y.fname = 'bodyPitch';  jpdf.(jname).y.lname = 'body pitch';
jpdf.(jname).bin.(jpdf.(jname).x.fname) = bin.(jpdf.(jname).x.fname);
jpdf.(jname).bin.(jpdf.(jname).y.fname) = bin.(jpdf.(jname).y.fname);
jpdf.(jname).data = sum(cat(3,iut{:}),3);
jpdf.(jname).label = {'Joint Probability Distribution', ...
                    ['between ',jpdf.(jname).x.lname,' and ',jpdf.(jname).y.lname],...
                    'during rearing in the arena perimeter',...
                    'head distance from center <30cm, maze radius: 40cm'};
figure();
imagesc(jpdf.(jname).bin.(jpdf.(jname).x.fname).centers,    ...
        jpdf.(jname).bin.(jpdf.(jname).y.fname).centers,    ...
        jpdf.(jname).data'                                  ...
        );
xlabel(bin.(jpdf.(jname).x.fname).label);
ylabel(bin.(jpdf.(jname).y.fname).label);
title(jpdf.(jname).label);
axis(gca(),'xy');
colormap(gca(),'jet');
set(gca(),'ColorScale','log');
cax = colorbar();
ylabel(cax,'Count');
saveas(gcf(),fullfile('/storage/gravio/figures/analysis/afig/',['jpdf_',jname,'_REAR_interior.pdf']));