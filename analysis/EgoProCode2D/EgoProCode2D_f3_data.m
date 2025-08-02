


hbaBin.edges = [-1.2,-0.2,0.2,1.2];
hbaBin.centers = mean([hbaBin.edges(1:end-1);hbaBin.edges(2:end)]);
hbaBin.count = numel(hbaBin.centers);
hbaBin.color = [0,1,0;...
                0,0,1;...
                1,0,0];
hbaBin.key = 'LCR';
hbaBin.label = {'Left','Center','Right'};


phzBin.edges = linspace(0.5,2*pi-0.5,4);
phzBin.centers = (phzBin.edges(1:end-1)+phzBin.edges(2:end))./2;
phzBin.count = numel(phzBin.centers);
phzBin.color = cool(phzBin.count);
phzBin.key = 'DTA';
phzBin.label = {'Descending','Trough','Ascending'};


hvfBin.edges = [-25,-5,5,25,80];
hvfBin.centers = (hvfBin.edges(1:end-1)+hvfBin.edges(2:end))./2;
hvfBin.count = numel(hvfBin.centers);
%hvfBin.color = summer(hvfBin.count);
hvfBin.color = lines(hvfBin.count);
hvfBin.key = 'RISF';
hvfBin.label = {'Reverse','Immobile','Slow','Fast'};
hvfBin.unit = 'cm/s';



EgoProCode2D_f2_data_egoHvf();
EgoProCode2D_f2_data_egoHvfPhz();


