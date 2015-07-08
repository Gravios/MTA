function stss = walk_ang(Trial,varargin)
%function sind = walk_lang(Trial)
%Segments the walking peroids by head angle
% 
%

[angThresh,mode] = DefaultArgs(varargin,{-.45,'threshold'});

stss = {};

switch mode
  case 'threshold'
    wind = Trial.stc{'w'}.copy;
    lang = Trial.stc{'w'}.copy;
    ang = Trial.ang.copy;
    ang.create(Trial,Trial.load('xyz'));
    lang.data = ThreshCross(ang(:,5,7,2)<angThresh,.5,20);
    lang.label = 'lang';
    lang.key = 'p';
    lang.filename = [Trial.filebase '.sst.lang.p.mat'];
    stss{1} = lang&wind.data;

    hang = Trial.stc{'w'}.copy;
    hang.data = ThreshCross(ang(:,5,7,2)>angThresh,.5,round(.25/ang.sampleRate));
    hang.label = 'hang';
    hang.key = 'c';
    hang.filename = [Trial.filebase '.sst.hang.c.mat'];
    stss{2} = hang&wind.data;

  case 'hmm'
    
end
