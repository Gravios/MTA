
MjgER2016_load_data();
configure_default_args();

phzCorrection = [pi*1.25,pi*1.25,                             ... er01
                 pi/2,pi/2,pi/2,                    ... ER06
                 pi/1.25,pi/1.25,                             ... Ed10
                 0,0,0,0,0,0,0,0,0,                 ... jg04                 
                 pi/4,pi/4,pi/4,pi/4,pi/4,pi/4,pi/4 ... jg05
                 pi/4,pi/4,pi,pi/1.25,pi*1.25]; % new units - jg05, jg05, ER06, Ed10, er01



pft = cf(@(t,u)  pfs_2d_theta(t,u,'overwrite',true),          Trials, units);
pfs = cf(@(t,u)  pfs_2d_states(t,u,'overwrite',true),         Trials, units);
bfs = cf(@(t,u)  compute_bhv_ratemaps(t,u,'overwrite',true),  Trials, units);
pftt = cf(@(t,u,tr,pc)                                        ...
          pfs_2d_thetaT(t,u,tr,pc,'overwrite',true),          ...
          Trials, units, num2cell([sessionList.thetaRefGeneral]), num2cell(phzCorrection));


t = 19;
[mxr,mxp] = maxRate(pft{t},units{t});
[mxrt,mxpt] = maxRate(pftt{t},units{t});
figure;
for u = units{t},
    subplot(121);    
    plot(pft{t},u,1,'text');
    Lines(mxp(u==units{t},1),[],'k');
    Lines([],mxp(u==units{t},2),'k');    
    subplot(122);
    plot(pftt{t},u,1,'text');    
    Lines(mxpt(u==units{t},1),[],'k');
    Lines([],mxpt(u==units{t},2),'k');    
    title(num2str(u));
    waitforbuttonpress();
end



t = 19;
[mxr,mxp] = maxRate(pft{t},units{t});
[mxrt,mxpt] = maxRate(pftt{t},units{t});
figure;
for u = units{t},
    subplot(131);    
    plot(bfs{t},u,1,'text',[],false);
    subplot(132);
    plot(pft{t},u,1,'text');    
    Lines(mxpt(u==units{t},1),[],'k');
    Lines([],mxpt(u==units{t},2),'k');    
    title(num2str(u));
    subplot(133);
    plot(pftt{t},u,1,'text');    
    Lines(mxpt(u==units{t},1),[],'k');
    Lines([],mxpt(u==units{t},2),'k');    
    title(num2str(u));
    waitforbuttonpress();
end


% This script was made to examine the posibility of using only the theta-trough spikes for theta place
% field computation.

% This may mean new place field centers (PFCs) for each unit which may result in better decoding 
% BUT Decoding may not be determined by rate probability.
% ALTERNATIVE : use mode PFC of the active placecell set using the mean-shift procedure by Fukunaga and hostetler (1975)
% MODIFICATION : use aformentioned mode with a weighted position history which has "momentum" -> refine definition




%figure,plot(mxrt-mxr,vecnorm(mxp-mxpt,2,2),'.');XB