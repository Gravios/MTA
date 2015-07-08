

function figure_controls_gui(hfig,index,indArray,mode)
hfig = get(hfig,'Parent');
eval([mode '(hfig,index,indArray)']);
end

function kb(hfig,index,indArray)
B = waitforbuttonpress;
whatkey = get(hfig,'CurrentCharacter');
if ~B,return,end
switch double(whatkey)
    case double('i')
        index = input('Enter index #: ');
        if sum(index==indArray),
            return,
        else
            index = input('Invalid Selection\nEnter index #: ');
        end
    case double('n')
        index = indArray(find(index==indArray)+1);
    case double('p')
        index = indArray(find(index==indArray)-1);
    case double('q')
        index = -1;
end
set(hfig,'Name',num2str(index));
end

function forwardButton_Callback(hfig,index,indArray)
    index = indArray(find(index==indArray)+1);
    set(hfig,'Name',num2str(index));
end


function backwardButton_Callback(hfig,index,indArray)
    index = indArray(find(index==indArray)-1);
    set(hfig,'Name',num2str(index));
end

function exitButton_Callback(hfig,index,indArray)
    index = -1;
    set(hfig,'Name',num2str(index));
end


function popupmenu_Callback(hfig,index,indArray)
    index = -1;
    set(hfig,'Name',num2str(index));
end


function printButton_Callback(hfig,index,indArray)
    if ispc, 
        userdir= getenv('USERPROFILE'); 
    else
        userdir= getenv('HOME');
    end
    dest = fullfile(userdir,'figures');
    folder = ['figs' datestr(now,29)];
    if ~exist(fullfile(dest,folder),'dir'),mkdir(fullfile(dest,folder));end
    filenamepng = fullfile(dest,folder,['f' datestr(now,30) '.png']);
    filenameeps = fullfile(dest,folder,['f' datestr(now,30) '.eps']);
    parent = hfig;%get(hfig,'Parent');
    print(parent,'-dpng',filenamepng);
    print(parent,'-dpsc2',filenameeps);
end
