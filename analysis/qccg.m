function qccg(Session,varargin)
[bhvs,filters] = DefaultArgs(varargin,{{'rear','walk'},{}});


bhv_on  = Session.Bhv.getState(bhvs{s}).state(:,1);
bhv_off = Session.Bhv.getState(bhvs{s}).state(:,2);

%[bhv_on,bhv_off] = Session.Bhv.exclude({'rear','bturn'},'walk',3);

for s = 1:numel(bhvs),
    MTAccg(Session,...
           bhvs{s},
           ['Testing a revision of MTAccg: added c_c_m feature'],...
           {bhv_on,bhv_off},...
           {[bhvs{s} '_onset'],[bhvs{s} '_offset']},
           [],...         Units to be processed 
           true,...       overwrite previously calculated values 
           false,...      add a six digit random number to the filename
           'abs_dist',... method used during partitioning
           4,...          number of partitions
           pfmp,...       partition feature used to distribute events into partitions
           1,...          keep a record of the xy position of the events
           0,...          surrogate sample size when selcting from a periods
           1);%           number of iterations
end