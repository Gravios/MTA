function [x,y] = estimateMazeCenter(Session)


switch Session.maze.shape
  case 'circle'

    xyz = Trial.load('xyz');
    xy = sq(xyz(:,'head_front',[1,2]));
    xy = xy(nniz(xy),:);
    mxy = mean(xy);
    xy = bsxfun(@minus,xy,mxy);

    initCo = [0,max(abs(diff(Session.maze.boundaries([1:2],:),1,2)))];
    rotStep = .1;
    sampSize = round(5*xyz.sampleRate);
    rotSteps = -pi:rotStep:pi;
    sind = [];
    pco = [];
    for i = 1:numel(rotSteps),
        rotCo = initCo*[cos(rotSteps(i)),-sin(rotSteps(i));sin(rotSteps(i)),cos(rotSteps(i))];
        [~,xyInd] = sort(sqrt(sum(bsxfun(@minus,xy,rotCo).^2,2)));
        sind(end+1:end+sampSize) = xyInd(1:sampSize);
        co = polyfit(xy(sind(end-sampSize:end),1),xy(sind(end-sampSize:end),2),1);
        pco(i,:) = [-co(1)^-1,(co(1)+co(1)^-1)+co(2)];
    end

    eCenter = []
    for i = 1:size(pco,1),
        for j = 1:size(pco,1),
            if i~=j,
                x = -(pco(i,2)-pco(j,2))/(pco(i,1)-pco(j,1));
                eCenter(end+1,:) = [x,polyval(pco(i,:),x)];
            end
        end
    end
    
                                
         
  case 'rectangle'
  case 'square'
  otherwise
    error('MTA:utilities:estimateMazeCenter:InvalidShape');
end
