
addpath(genpath('/gpfs01/sirota/bach/homes/gravio/modules/fieldtrip-20140510/'));
cd /gpfs01/sirota/bach/data/gravio/nlx/jg05-20120312-nlx/

cfg         = [];
cfg.dataset = 'jg05-20120312-01-cof-nlx-065.ncs';
data65        = ft_preprocessing(cfg);


ts1 = ft_read_data('jg05-20120312-01-cof-nlx-001.ncs','timestamp',true);
ts2 = ft_read_data('jg05-20120312-01-cof-nlx-002.ncs','timestamp',true);
ts65 = ft_read_data('jg05-20120312-01-cof-nlx-065.ncs','timestamp',true);
ts90 = ft_read_data('jg05-20120312-01-cof-nlx-090.ncs','timestamp',true);

fname = {'jg05-20120312-01-cof-nlx-001.ncs','jg05-20120312-01-cof-nlx-002.ncs'};
[data_2chan] = ft_read_neuralynx_interp(fname);


ts1 = double(ts1);
ts2 = double(ts2);
ts65 = double(ts65);
ts90 = double(ts90);

(ts1d(end)-ts1d(1))/10^6/60/60



ts1d = diff(ts1);
ts1g = ts1d(ts1d>31);
numel(ts1g)

ts2d = diff(ts2);
ts2g = ts1d(ts2d>31);
numel(ts2g)


ts90d = diff(ts90);
ts90g = ts90d(ts90d>31);
numel(ts90g)


xs = cell(1,2);
[xs{:}] = meshgrid(1:1000,1:1000);
xs = reshape(cat(3,xs{:}),[],2);
%xs = cat(3,xs{:});
cs = repmat([30,31],size(xs,1),1);
mm = sum(xs.*cs,2);


ts90g(10)

astat = zeros([1,10000]);
icount = 1;
for i = round(linspace(1,size(ts90d,2)-11000,10000)),
astat(icount) = sum(ts90d([1+i]:[1000+i])==30)/sum(ts90d([1+i]:[1000+i])==31);
icount=icount+1;
end




a = mean(astat);

i =1:10;
y = round(ts90g(i)/(31+30*a))
x = (ts90g(i)-y*31)/30

xs(find(mm==ts90g(3)),:)*[30;31]

[x,y]

gapIndSize = 

xs(find(mm==ts90g(3)),:)

numel(find(mm==ts90g(10)))