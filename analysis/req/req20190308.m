


% HISTOGRAM : 2D : JPDF : forward versus lateral relative to head position
%             subdivided by Trial and State
%             get smoothed center either: 
%                 1. convolved local mininima
%                 2. gaussian fit
%
% HISTOGRAM : 2D : JPDF : forward versus lateral relative to head position
%             all Trials subdivided by state and angle between head trajectory vector 
%             and head direction vector
%             
%

% dcer: xy error - distance between posterior maximum and current position
% 


MjgER2016_load_data();

sampleRate = 30;   % Hz [C:250]
spikeWindow = 0.3; % ms [C:0.03]
states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};

sids = 17:23;

[posEstComXYH,posEstMaxXYH,posEstSaxXYH,posteriorMaxXYH] = ...
    cf(@(t,u)  bhv_decode(t,sampleRate,[],u,'xyh',[],[],spikeWindow,false), Trials(sids),units(sids));
[posEstComXYB,posEstMaxXYB,posEstSaxXYB,posteriorMaxXYB] = ...
    cf(@(t,u)  bhv_decode(t,sampleRate,[],u,'xyb',[],[],spikeWindow,false), Trials(sids),units(sids));

pMax = posteriorMaxXYH;
stc  = cf(@(t)    t.load('stc','msnn_ppsvd_raux'),               Trials(sids)   );
xyz  = cf(@(t)    resample(preproc_xyz(t,'trb'),sampleRate),     Trials(sids)   );
fxyz = cf(@(x)    filter(copy(x),'ButFilter',3,10,'low'),        xyz            );
fet  = cf(@(t)    fet_HB_pitchB(t,sampleRate),                   Trials(sids)   );
stcm = cf(@(s,x)  stc2mat(s,x,states),                           stc, xyz       );

tvec = cf(@(x)    circshift(x(:,'hcom',[1,2]),...
                            -round(sampleRate.*0.05))...
                  -x(:,'hcom',[1,2]),                            fxyz           );
tvec = cf(@(x)    sq(bsxfun(@rdivide,x,sqrt(sum(x.^2,3)))),      tvec           );
tvec = cf(@(x)    cat(3,x,sq(x)*[0,-1;1,0]),                     tvec           );

hvec = cf(@(x)    x(:,'head_front',[1,2])-x(:,'head_back',[1,2]),fxyz           );
hvec = cf(@(x)    sq(bsxfun(@rdivide,x,sqrt(sum(x.^2,3)))),      hvec           );
hvec = cf(@(x)    cat(3,x,sq(x)*[0,-1;1,0]),                     hvec           );

anghd = circ_dist(atan2(t(:,1),pfmr(:,2)),atan2(pfhr(:,1),pfhr(:,2)));

%drz  = cf(@(t))

% COUNT of active units within spiking window
uInc = cf(@(t,u,x)                                                                               ...
          sum(get(load(t,'ufr',x,[],u,spikeWindow,true,'gauss'),'data')>0.2,2),                  ...
          Trials(sids), units(sids), xyz);

% ERROR of bayesian decoding within head reference frame
dcer = cf(@(peXYH,peXYB,hv,x,f)                                                                  ...
          [multiprod(peXYH(:,[1,2])-sq(x(:,'hcom',[1,2])),hv,2,[2,3]),                           ...
          f(:,1)-peXYH(:,3),f(:,2)-peXYB(:,3)],                                                  ...
          posEstSaxXYH,posEstSaxXYB,hvec,xyz,fet);



ind = cf(@(ui,pm,st)  ui>=4 & pm>0.01 & st(:,1) & any(st(:,[2,3,4,5,6]),2),  uInc, pMax, stcm);
% HISTOGRAM : 2D : JPDF : forward versus lateral relative to head position ------------------------


% LOWRES 
edx = linspace(-900,900,100);
edy = linspace(-900,900,100);
edh = linspace(-2,2,50);
edb = linspace(-2,2,50);
% HIGHRES 
edx = linspace(-300,300,50);
edy = linspace(-300,300,50);
edh = linspace(-2,2,50);
edb = linspace(-2,2,50);

pageWidth  = 21.0;
pageHeight = 29.7;

hfig = figure(666005);
clf(hfig);
hfig.Units = 'centimeters';
hfig.Position = [1, 1, pageWidth,pageHeight];
hfig.PaperPositionMode = 'auto';

sp = tight_subplot(numel(sids)+1,6,0.01,0.1,0.1);
sp = reshape(sp,[6,numel(sids)+1]);
minds = {};
for t = 1:numel(sids)+1,
    y = t;    
    if t == numel(sids)+1,
        t = 1:numel(sids);
    end        
    for s = 1:6,
        x = s;
        ind = cf(@(ui,pm,st)  ui>=4 & pm>0.01 & st(:,1) & any(st(:,[s]),2),  uInc(t), pMax(t), stcm(t));
        jpdfXY = cf(@(e,i)  hist2(e(i,[1,2]),edx,edy), dcer(t), ind);
        %jpdfHB = cf(@(e,i)  hist2(e(i,[3,4]),edh,edb), dcer(t), ind);
        axes(sp(x,y));
        imagesc(edx,edy,sum(cat(3,jpdfXY{:}),3)'./sum(reshape(cat(3,jpdfXY{:}),[],1)));
        %minds{y,s} = LocalMinima2(-sum(cat(3,jpdfXY{:}),3),0,100);
        axis('xy');
        Lines([],0,'r',1);
        Lines(0,[],'r',1);
        if y==1,  title(states{s}); end
        if x==1&&y~=numel(sids)+1,
            ylabel(Trials{sids(t)}.name); 
        elseif x==1&&y==numel(sids)+1,
            ylabel('jg05'); 
        end        
% $$$ subplot(122);
% $$$ imagesc(edh,edb,sum(cat(3,jpdfHB{:}),3)');
% $$$ axis('xy');
        set(gca,'XTickLabels',{});
        set(gca,'YTickLabels',{});
    end
end
ForAllSubplots('caxis([0,0.01])');
af(@(s) set(s,'Units','centimeters'), sp);
af(@(s) set(s,'Position',[s.Position(1:2),2.5,2.5]), sp);
drawnow();
pause(0.3);
cb = colorbar(sp(end));
cb.Position(1) = cb.Position(1)+0.075;
ylabel(cb,'Probability')

print(gcf,'-depsc2',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['pp_pop_drz_hrz.eps']]);
print(gcf,'-dpng',...
      ['/storage/share/Projects/BehaviorPlaceCode/phase_precession/',...
       ['pp_pop_drz_hrz.png']]);


hfig.Units = 'centimeters';


[x,y] = cf(@(m) ind2sub([99,99],m), mins)


% DIAGNOSTIC : check real and decoded trajectories 

figure();
s = 4;
sp = gobjects([0,1]);
sp(end+1) = subplot(4,1,1);
hold('on');
plot(xyz{s}(ind{s},'hcom',1));
plot(posEstSaxXYH{s}(ind{s},1),'r');
plot(posEstSaxXYB{s}(ind{s},1),'m');
sp(end+1) = subplot(4,1,2);
hold('on');
plot(xyz{s}(ind{s},'hcom',2));
plot(posEstSaxXYH{s}(ind{s},2),'r');
plot(posEstSaxXYB{s}(ind{s},2),'m');
sp(end+1) = subplot(4,1,3);
hold('on');
plot(fet{s}(ind{s},1));
plot(posEstSaxXYH{s}(ind{s},3),'r');
sp(end+1) = subplot(4,1,4);
hold('on');
plot(fet{s}(ind{s},2));
plot(posEstSaxXYB{s}(ind{s},3),'m');
linkaxes(sp,'x');

