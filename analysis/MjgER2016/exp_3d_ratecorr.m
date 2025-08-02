

Trial = MTATrial.validate('ER06-20130614.cof.all');
Stc = Trial.load('stc','default');

units = select_units(Trial,18,'pyr');
overwrite = false;

Trial = MTATrial.validate('ER06-20130614.cof.gnd');
Trial.maze.boundaries(end) = 450;
Trial.stc = Stc;
Trial.load('stc');

pfs{1} = MTAApfs(Trial,units,'theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);



Trial = MTATrial.validate('ER06-20130614.cof.fly');
Trial.maze.boundaries(end) = 450;
Trial.stc = Stc.copy;
Trial.load('stc');

if isempty(Trial.stc.gsi('f')),Trial = labelFlight(Trial);end
fstc = Trial.load('stc','flight');

pZ(Trial)
xyz = Trial.load('xyz');
hold on,plot(xyz(:,1,3));
Lines(fstc{'f'}(:),[],'r');



pfs{2} = MTAApfs(Trial,units,'theta-flight',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);

pfs{3} = MTAApfs(Trial,units,'flight&theta',...
                 overwrite,...
                 'numIter',1,...
                 'binDims',[20,20,20],...
                 'type','xyz',...
                 'SmoothingWeights',[2.2,2.2,2.2]);

