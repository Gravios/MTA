
%% Compare rearss with rfccg

figure(102)
set(gcf,'CurrentCharacter','l');
unit = 1;
while 1,
    clf
    subplot(211)
    bar(tbin,rfccg(:,2,unit));axis tight
    subplot(212)
    bar(tbin,ConArgs{1}(1,:,2,unit));axis tight
    waitforbuttonpress
    whatkey = get(gcf,'CurrentCharacter');
    switch double(whatkey)
      case double('i')
        unit = input('Enter unit #: ');
      case double('n')
        unit = unit+1;
      case double('p')
        unit=unit-1;
      case double('q')
        return
    end
end
