% req20171216 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: parameter scripts for testing locate_rhm_ncp_coherence_point
%  Bugs: NA

FigPath = '/storage/gravio/figures/analysis/locate_max_rhmXncp_coherence';
create_directory(FigPath);

sessionList = get_session_list('ncp');

% INITIAL search space
%    i = [-50:5:50];   longitudinal
%    j = [-50:5:50];   lateral
%    k = [-50:5:50]; vertical 

% Ed05-20140529.ont.all
cvals = [min(vxyz(:)),max(vxyz(:))];

hfig = figure();
hfig.Position = [1,50,1000,700];
hfig.PaperPositionMode = 'auto';
offset = 16;

for z = 1:5,
    for m = 1:7,
        hax = subplot2(5,7,z,m);
        imagesc(i,j,sq(vxyz(m,:,:,z+offset))');
        caxis(cvals);
        daspect([1,1,1])

        title(['Z disp: ',num2str(k(z+offset)),'mm']);                
        if m == 1,
            ylabel('Y disp (mm)');
            title({nhm{m+1},['Z disp: ',num2str(k(z+offset)),'mm']});            
        end

        if z == 1,
            title({nhm{m+1},['Z disp: ',num2str(k(z+offset)),'mm']});            
        end

        if z == 5,
            xlabel('X disp (mm)');
        end

        if z == 5 && m == 7,
            cax = colorbar();
            cax.Position([1,2,4]) = [cax.Position(1)+0.04,hax.Position([2,4])];
        end

    end
end
suptitle(['mean rhm ncp cohorence for displaced hcoms: ',Session.filebase])

% SAVE figure
% SAVE png
FigName = 'prot_Ed05-20140529.ont.all_search_headRefSpace_max_rhmXncp_coherence';
savefig(hfig,fullfile(FigPath,[FigName,'.fig']))
print(hfig,'-depsc2',fullfile(FigPath,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(FigPath,[FigName,'.png']));



hfig = figure();
plot([-150:5:150],vz);
xlabel('Z disp from trb hcom (mm)');
ylabel('rhmXncp Coherence');
title('Ed20140529.ont.all');
FigName = 'prot_Ed05-20140529.ont.all_search_headRefSpace_Z_max_rhmXncp_coherence';
savefig(hfig,fullfile(FigPath,[FigName,'.fig']))
print(hfig,'-depsc2',fullfile(FigPath,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(FigPath,[FigName,'.png']));



% ADJUST search space
%    i = [-20:10:20];   longitudinal
%    j = [-20:10:20];   lateral
%    k = [-150:10:150]; vertical 

Session = MTASession.validate(sessionList(5)); % 'Ed05-20140529.ont.all'
Session = MTASession.validate(sessionList(4)); % 'Ed05-20140529.ont.all'
Session = MTASession.validate(sessionList(2)); % 'Ed03-20140709.cof.all'


cvals = [min(vxyz(:)),max(vxyz(:))];

hfig = figure();
hfig.Position = [1,50,1000,700];
hfig.PaperPositionMode = 'auto';
offset = 35;

for z = 1:5,
    for m = 1:7,
        hax = subplot2(5,7,z,m);
        imagesc(i,j,sq(vxyz(m,:,:,z+offset))');
        caxis(cvals);
        daspect([1,1,1])

        title(['Z disp: ',num2str(k(z+offset)),'mm']);                
        if m == 1,
            ylabel('Y disp (mm)');
            title({nhm{m+1},['Z disp: ',num2str(k(z+offset)),'mm']});            
        end

        if z == 1,
            title({nhm{m+1},['Z disp: ',num2str(k(z+offset)),'mm']});            
        end

        if z == 5,
            xlabel('X disp (mm)');
        end

        if z == 5 && m == 7,
            cax = colorbar();
            cax.Position([1,2,4]) = [cax.Position(1)+0.04,hax.Position([2,4])];
        end

    end
end
suptitle(['mean rhm ncp cohorence for displaced hcoms: ',Session.filebase])

% SAVE figure
% SAVE png
FigName = [Session.filebase,'_search_headRefSpace_max_rhmXncp_coherence'];
savefig(hfig,fullfile(FigPath,[FigName,'.fig']))
print(hfig,'-depsc2',fullfile(FigPath,[FigName,'.eps']));
print(hfig,'-dpng',  fullfile(FigPath,[FigName,'.png']));

