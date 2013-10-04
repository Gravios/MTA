%ds = load ('/data/homes/gravio/data/analysis/jg05-20120310/jg05-20120310.cof.all.frameset_12505_14235.mat')
ds = load ('/data/homes/gravio/data/analysis/jg05-20120310/jg05-20120310.cof.all.frameset_12498_14223.mat')
%ds = load ('/data/homes/gravio/data/analysis/jg05-20120309/jg05-20120309.cof.all.frameset_50331_51147.mat')
frameset = {};
for i=1:length(ds.record),
    frameset {i} = rmfield(ds.record{i},'index');
end

tmovie = frameset{1};
for i=1:length(frameset),
tmovie(i) =  frameset{i};
tmovie(i).cdata = tmovie(i).cdata(81:650,101:800,:);
end


movie2avi(tmovie,'ant_rat_bhv2.avi','compression','None')

%% Avifile Test
aviobj = avifile('ant_avi_test.avi');
for i = 1:length(frameset),
    aviobj=addframe(aviobj,frameset{i});
end,
aviobj = close(aviobj);