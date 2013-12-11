

function figure_controls(hfig,index,indArray,mode)
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