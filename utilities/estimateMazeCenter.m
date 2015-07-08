function [x,y] = estimateMazeCenter(Session)


switch Session.maze.shape
  case 'circle'

    xyz = Session.load('xyz');
    xy = sq(xyz(:,'head_front',[1,2]));
    xy = xy(nniz(xy),:);
    mxy = mean(xy);
    xy = bsxfun(@minus,xy,mxy);
    
    radius = 420;

    initCo = [0,max(abs(diff(Session.maze.boundaries([1:2],:),1,2)))];
    rotStep = .1/pi;
    sampSize = round(10*xyz.sampleRate);
    rotSteps = -pi:rotStep:3*pi;

    sind = [];
    pco = [];
    
    for i = 1:numel(rotSteps),
        rotCo = initCo*[cos(rotSteps(i)),-sin(rotSteps(i));sin(rotSteps(i)),cos(rotSteps(i))];
        [~,xyInd] = sort(sqrt(sum(bsxfun(@minus,xy,rotCo).^2,2)));
        sind(end+1:end+sampSize) = xyInd(1:sampSize);
        co = polyfit(xy(sind(end-sampSize+1:end),1),xy(sind(end-sampSize+1:end),2),1);
        pco(i,:) = [-co(1)^-1,mean(xy(sind(end-sampSize+1:end),2))*2-(co(1)-co(1)^-1)*mean(xy(sind(end-sampSize+1:end),1))-co(2)];
    end

    eCenter = [];
    for i = 1:size(pco,1),
        for j = 1:size(pco,1),
            if i~=j,
                x = -(pco(i,2)-pco(j,2))/(pco(i,1)-pco(j,1));
                eCenter(end+1,:) = [x,polyval(pco(i,:),x)];
            end
        end
    end
    
    xyz.data = [xyz(:,:,1)-(median(eCenter(:,1))+mxy(1)),xyz(:,:,2)-(median(eCenter(:,2))+mxy(2)),xyz(:,:,3)];
    
    
    
    figure,plot(xyz(:,7,1)-(median(eCenter(:,1))+mxy(1)),xyz(:,7,2)-(median(eCenter(:,2))+mxy(2)),'.')
    ang = 0:0.01:2*pi;    
    x = radius*cos(ang);
    y = radius*sin(ang);
    hold on,plot(x,y,'r')

    oscore = sum((xyz(:,7,1).^2+xyz(:,7,2).^2 - radius.^2)>0);
    
    escore = 
    
    figure,plot(xyz(:,7,1),xyz(:,7,2),'.')
    hold on,plot(x,y,'r')
    figure,plot(eCenter(:,1)+mxy(:,1),eCenter(:,2)+mxy(:,2),'.')    
    
%     figure,plot(eCenter(:,1),eCenter(:,2),'.')    
%     hold on
%     plot(xy(unique(sind),1),xy(unique(sind),2),'r.')
%     plot(xy(sind(end-sampSize+1:end),1),xy(sind(end-sampSize+1:end),2),'g.')
%     plot(-1000:1000,polyval(pco(end,:),-1000:1000),'g')
% 
%     plot(xy(sind(sampSize*30+1:sampSize*31),1),xy(sind(sampSize*30+1:sampSize*31),2),'m.')
%     plot(-1000:1000,polyval(pco(30,:),-1000:1000),'m')
%     for i = 1:2:20,
%         plot(-1000:1000,polyval(pco(i,:),-1000:1000))
%     end
%     
    
         
  case 'rectangle'
  case 'square'
  otherwise
    error('MTA:utilities:estimateMazeCenter:InvalidShape');
end
