function ACh = load_all(Trial),

signalType = 'ACh';
fbrLoc = regexp(af(@(f)  f.name,  dir(fullfile(Trial.spath,['*.fbr.*.',signalType]))),...
           ['(?<=fbr\.)\w+(?=\.',signalType,'$)'],'once','match')';
ACh = cf(@(loc) Session.load('fbr', loc, 'ACh', 'ACh'), fbrLoc);
nt  = cf(@(loc) Session.load('fbr', loc, 'NULL', 'ACh'), fbrLoc);
for a = 2:numel(ACh),  ACh{a}.data(:,end+1) = nt{a}.data;  end
for a = 2:numel(ACh),  ACh{1}.data(:,:,end+1) = ACh{a}.data;  end
ACh = ACh{1};
