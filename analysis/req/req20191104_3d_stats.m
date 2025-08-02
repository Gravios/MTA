



ctrl = load(fullfile(MTA_PROJECT_PATH,'analysis',...
                     ['MjgER2016_req20191104_jpdf_hbaCorrected_3d_19cb8a5111c2a4aeac651449bf491daa.mat']));

perm = load(fullfile(MTA_PROJECT_PATH,'analysis',...
                     ['MjgER2016_req20191104_jpdf_hbaCorrected_3d_perm19cb8a5111c2a4aeac651449bf491daa.mat']));




g = 1;
h = 5;
d = 6;
k = 0;
for s = [2,3,4];
    %for s = [1,6,7,8];
    k = k+1;        
    hfig = figure();
    hfig.Units = 'centimeters';
    hfig.Position = [(k-1).*18,5,18,24];
    
    for z = 1:3,
        for p = 1:8,
            for e = 1:2
                subplot2(8,nbins*2,p,z+(e-1).*3);
                nmask = sq(double(ctrl.sCnt(z,:,:,2,g,h,d,s)>10000));
                nmask(~nmask) = nan;
                if e==1,
                    ca = [-60,120];
                else
                    ca = [-60, 60];
                end
% $$$         imagescnan({vbc{grps{g}(1)},vbc{grps{g}(2)},sq(smJpdf(2,phzOrder(p),:,:,e,g,s,:))'.*nmask'},...
% $$$                    ca,'linear',false,'colorMap',@jet);
% $$$         axis('xy');
                imagesc(ctrl.vbc{h},ctrl.vbc{d},sq(ctrl.smJpdf(2,phzOrder(p),z,:,:,e,g,h,d,s,:))'.*nmask');
                colormap('jet');
                caxis(ca)
                axis('xy');
                if e==1&&p==1&&z==1,
                    title(stlbls{s});
                end
            end
        end
    end
end
%%%>>>
