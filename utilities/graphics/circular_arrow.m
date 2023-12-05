function [hplots, arrows] = circular_arrow(axesHandle, radius, center, arrow_angle,angle, ...
                             direction, color, head_size, head_style)


xc = center(1);
yc = center(2);

xt = xc + radius;
yt = yc;

x1 = (xt-xc) * cos( arrow_angle + angle/2) - ...
     (yt-yc) * sin( arrow_angle + angle/2) + xc;
x2 = (xt-xc) * cos( arrow_angle - angle/2) - ...
     (yt-yc) * sin( arrow_angle - angle/2) + xc;
x0 = (xt-xc) * cos( arrow_angle) - ...
     (yt-yc) * sin( arrow_angle) + xc;

y1 = (xt-xc) * sin( arrow_angle + angle/2) - ...
     (yt-yc) * cos( arrow_angle + angle/2) + yc;
y2 = (xt-xc) * sin( arrow_angle - angle/2) - ...
     (yt-yc) * cos( arrow_angle - angle/2) + yc;
y0 = (xt-xc) * sin( arrow_angle) - ...
     (yt-yc) * cos( arrow_angle) + yc;

i = 1;

P1 = struct([]);
P2 = struct([]);

P1{1} = [x1;y1];
P1{2} = [x2;y2];

P2{1} = [x0;y0];
P2{2} = [x0;y0];

center = [xc; yc];
n = 1000;
v = struct([]);

hplots = gobjects([2,0]);

for i = 1:numel(P1)
    v1 = P1{i} - center;
    v2 = P2{i} - center;
    c = det([v1,v2]);
    a = linspace(0, atan2( abs(c), dot(v1,v2)), n);
    v3 = [0, -c; c, 0] * v1;
    v{i} = v1 * cos(a) + ((norm(v1)/norm(v3)) * v3) * sin(a);
    hplots(1,i) = plot(v{i}(1,:) + xc, ...
                       v{i}(2,:) + yc, ...
                       'Color', color);
end

position = struct([]);
if direction == 1
    position{1} = [x2, y2, x2 - (v{2}(1,2)+xc), y2-(v{2}(2,2)+yc)];
elseif direction == -1
    position{1} = [x1, y1, x1 - (v{1}(1,2)+xc), y1-(v{1}(2,2)+yc)];
elseif direction == 2
    position{1} = [x2, y2, x2 - (v{2}(1,2)+xc), y2-(v{2}(2,2)+yc)];
    position{2} = [x1, y1, x1 - (v{1}(1,2)+xc), y1-(v{1}(2,2)+yc)];
end


arrows = gobjects([0,1]);
% arrow heads
for i = 1:numel(position)
    arrows(i) = annotation('arrow');
    set(  arrows(i),               ...
           'Parent', axesHandle,   ...
         'Position', position{i},  ...
       'HeadLength', head_size(2), ...
        'HeadWidth', head_size(1), ...
        'HeadStyle', head_style,   ...
        'LineStyle', 'none',       ...
            'Color', color         ...
        );
end
