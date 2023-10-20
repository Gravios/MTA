

Trial.stc.states = {};
% LABEL Non nan periods
gper = ThreshCross(double(rat(:,1,1)~=eps),0.5,10);

Trial.stc.addState(Trial.spath,                     ... path to project folder
                   Trial.filebase,                  ... filebase
                   gper,                            ... state periods [n,2] 
                   rat.sampleRate,                  ... period sample rate
                   Trial.sync.copy,                 ... Trial synchronization object
                   Trial.sync.data(1),              ... Trial synchronization origin
                   'gper',                          ... state label
                   'a');
Trial.stc.states{end}.save(true);

% LABEL speed thresholded periods VEL
ratFilt = filter(copy(ratAC),'ButFilter',4,1,'low');
%figure();plot(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10)
vper = ThreshCross(double(sqrt(sum(diff(ratFilt(:,1,1:2)).^2,3)).*ratFilt.sampleRate./10>2),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   vper,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'vel','v');
Trial.stc.states{end}.save(true);


% LABEL Height thresholded periods REAR
ratFiltLow = filter(copy(ratAC),'ButFilter',4,0.75,'low');
rper = ThreshCross(double(ratFiltLow(:,1,3)>150),0.5,10);
Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   rper,...
                   rat.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   'rear','r');
rper = Trial.stc.states{end}.copy();
rper = rper+[-0.2,0.2];
rper.clean();
rper.save(true);

% $$$ Trial.stc.addState(Trial.spath,                     ... path to project folder
% $$$                    Trial.filebase,                  ... filebase
% $$$                    [bTimes-10*rat.sampleRate,bTimes],                            ... state periods [n,2] 
% $$$                    rat.sampleRate,                  ... period sample rate
% $$$                    Trial.sync.copy,                 ... Trial synchronization object
% $$$                    Trial.sync.data(1),              ... Trial synchronization origin
% $$$                    'approach',                          ... state label
% $$$                    'y');
% $$$ Trial.stc.states{end}.save(1);
% $$$ 
% $$$ 
% $$$ Trial.stc.addState(Trial.spath,                     ... path to project folder
% $$$                    Trial.filebase,                  ... filebase
% $$$                    [bTimes,bTimes+10*rat.sampleRate],                            ... state periods [n,2] 
% $$$                    rat.sampleRate,                  ... period sample rate
% $$$                    Trial.sync.copy,                 ... Trial synchronization object
% $$$                    Trial.sync.data(1),              ... Trial synchronization origin
% $$$                    'depart',                          ... state label
% $$$                    'd');
% $$$ Trial.stc.states{end}.save(1);



% $$$ Trial.stc.addState(Trial.spath,                     ... path to project folder
% $$$                    Trial.filebase,                  ... filebase
% $$$                    [bTimes-10*rat.sampleRate,bTimes],                            ... state periods [n,2] 
% $$$                    rat.sampleRate,                  ... period sample rate
% $$$                    Trial.sync.copy,                 ... Trial synchronization object
% $$$                    Trial.sync.data(1),              ... Trial synchronization origin
% $$$                    'approach',                          ... state label
% $$$                    'y');
% $$$ Trial.stc.states{end}.save(1);
% $$$ 
% $$$ 
% $$$ Trial.stc.addState(Trial.spath,                     ... path to project folder
% $$$                    Trial.filebase,                  ... filebase
% $$$                    [bTimes,bTimes+10*rat.sampleRate],                            ... state periods [n,2] 
% $$$                    rat.sampleRate,                  ... period sample rate
% $$$                    Trial.sync.copy,                 ... Trial synchronization object
% $$$                    Trial.sync.data(1),              ... Trial synchronization origin
% $$$                    'depart',                          ... state label
% $$$                    'd');
% $$$ Trial.stc.states{end}.save(1);
