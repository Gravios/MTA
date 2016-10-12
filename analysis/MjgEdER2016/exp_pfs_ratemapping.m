
Trial = MTATrial.validate('Ed10-20140817.cof.all');
Stc = Trial.load('stc','NN0317R');

Trial = 'Ed10-20140817.cof.gnd';

Trial = MTATrial.validate(Trial);
OwnDir = '/storage/gravio/ownCloud/MjgEdER2016/';
FigDir = 'exp_pfs_ratemapping';
mkdir(fullfile(OwnDir,FigDir))


pft = pfs_2d_theta(Trial,[],[],true);
mrt = pft.maxRate;
% Reduce clu list based on theta pfs max rate
%units = pft.data.clu(mrt>1);

%unitt = select_units(Trial,18);
Trial.load('nq');
% $$$ units = units(Trial.nq.SNR(units)>.75);
units = find(Trial.nq.AmpSym&mrt>1)';
