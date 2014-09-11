ses = MTASession([]);
ds = dir(ses.path.data);
snames = {ds(cellfun(@numel,{ds.name})==13).name};
for s =1:numel(snames),
    disp(snames{s})
try
    ts = dir(fullfile(ses.path.data,snames{s}));
    tnames = ts(~cellfun(@isempty,regexpi({ts.name},'trl'))).name;
    disp(tnames)
end
end
