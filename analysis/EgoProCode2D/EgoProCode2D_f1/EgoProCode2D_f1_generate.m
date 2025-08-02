

global MTA_PROJECT_PATH
partsPath = fullfile(fullfile(MTA_PROJECT_PATH,'analysis','EgoProCode2D','EgoProCode2D_figure_parts'));
overwrite = false;
%EgoProCode2D_f1_parts( overwrite );

% xmpl,plcfld, egofld allofild

% $$$ document = parse_inkscape_layout(fullfile(MTA_PATH,'analysis','EgoProCode2D',['EgoProCode2D_figure1.layout']),'EgoProCode2D_f1');
% $$$ document = parse_inkscape_layout(fullfile(MTA_PATH,'analysis','EgoProCode2D',['EgoProCode2D_figure1_alt1.layout']),'EgoProCode2D_f1');
document = parse_inkscape_layout(fullfile(MTA_PATH,'analysis','EgoProCode2D',['EgoProCode2D_figure1_alt2.layout']),'EgoProCode2D_f1');

document.body
document.subplots


%%%<<< Generate Figure
[hfig,fig,fax,sax] = set_figure_layout(figure(666001),'A4','portrait',[],2,2,1,0.6);
xlim(fax,[0,1]);ylim(fax,[0,1]);
hfig.Units = 'normalized'
widthScale  = document.body.width  ;
heightScale = document.body.height ;
k = 1;
for subplotName = fieldnames(document.subplots)'
    subplotName = subplotName{1};
    tempfig = openfig(fullfile(partsPath,[subplotName,'.fig']));
    tempfig.Units='normalized';
    ax = copyobj(tempfig.Children,hfig);
    % GET subplot handle
    sax = findobj(ax,'Tag',subplotName)
    % REPOSITION subplot
    for a = 1:numel(sax)
        sax(a).Position = [document.subplots.(subplotName).x/widthScale, ...
                        1 - document.subplots.(subplotName).y/heightScale - document.subplots.(subplotName).height/heightScale,...
                        document.subplots.(subplotName).width/widthScale,...
                        document.subplots.(subplotName).height/heightScale];
    end
    
    close(tempfig);
    
    disp(subplotName);
    k = k + 1;
end
%%%>>>



%%%<<< Render Individual Subplots
hfig = figure(666001);
render_subplot_by_name(hfig, ...
                       document, ...
                       partsPath, ...
                       'EgoProCode2D_f1_subplot_lateral_field_distrib_descending_phase')

subplotName = 'EgoProCode2D_f1_subplot_forward_field_distrib_every_phase'
subplotName = 'EgoProCode2D_f1_subplot_lateral_field_distrib_every_phase'
subplotName = 'EgoProCode2D_f1_ego_maxrate_stats_ca1';
subplotName = 'EgoProCode2D_f1_ego_maxrate_stats_ca3';
subplotName = 'EgoProCode2D_f1_ego_field_size_stats_ca1';
subplotName = 'EgoProCode2D_f1_ego_field_size_stats_ca3';
%subplotName = 'EgoProCode2D_f1_subplot_lateral_field_distrib_descending_phase';
hfig = figure(666001);
render_subplot_by_name(hfig, ...
                       document, ...
                       partsPath, ...
                       subplotName);


hfig = figure(666001);
render_subplot_by_name(hfig, ...
                       document, ...
                       partsPath, ...
                       'EgoProCode2D_f1_subplot_lateral_distrib_stats')

hfig = figure(666001);
render_subplot_by_name(hfig, ...
                       document, ...
                       partsPath, ...
'EgoProCode2D_f1_subplot_egofieldExample');
%'EgoProCode2D_f1_subplot_placefieldExample');
%%%>>>

