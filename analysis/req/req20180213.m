
Trial = Trials{19};
pft = pfs_2d_theta(Trial);
pfd = MjgER2016_drzfields(Trial,[],true);

units = select_placefields(Trial);

drz = compute_drz(Trial,units,pft);
xyz = preproc_xyz(Trial,'trb');
pch = fet_HB_pitchB(Trial,[],[],[],'Ed05-20140529.ont.all');
ang = create(MTADang,Trial,xyz);
stc = Trial.stc.copy();

spk = Trial.spk.copy();
spk.create(Trial,xyz.sampleRate,'theta',units,'deburst');    

nx = 1;
ny = numel(pfd)+1;
maxPfsRate = max(cell2mat(cf(@(p,u) p.maxRate(u,false),cat(2,{pft},pfd),repmat({units},[1,ny]))),[],2);

figure,

for u = units(:)'
    clf();
    mrate = maxPfsRate(u==units);
    subplot2(ny,nx,1,1);hold('on');
    plot(pft,u,'mean',true,mrate,false,0.99);
    plot(xyz(stc{'r'},'nose',1),xyz(stc{'r'},'nose',2),'.r');
    title(num2str(u));
    res = spk(u);
    for i = 1:ny-1
        subplot2(ny,nx,i+1,1);
        plot(pfd{i},u,'mean',true,mrate,false,0.99);
        hold('on');
        if     i==1, plot(drz(:,u==units),xyz(:,'nose',3),'.c','MarkerSize',1);
                     plot(drz(res,u==units),xyz(res,'nose',3),'.m','MarkerSize',5);
        elseif i==3, plot(drz(:,u==units),pch(:,1),'.c','MarkerSize',1);
                     plot(drz(res,u==units),pch(res,1),'.m','MarkerSize',1);                        
        elseif i==4, plot(drz(:,u==units),pch(:,2),'.c','MarkerSize',1);
                     plot(drz(res,u==units),pch(res,2),'.m','MarkerSize',5);
        end
    end
    drawnow();
    waitforbuttonpress();
end




figure,
plot(circ_dist(ang(:,'hcom','nose',2),ang(:,'spine_middle','spine_upper',2))+0.01);
hold('on');
plot(ang(:,'hcom','nose',2)-ang(:,'spine_middle','spine_upper',2));
plot(ang(:,'spine_middle','spine_upper',2));


