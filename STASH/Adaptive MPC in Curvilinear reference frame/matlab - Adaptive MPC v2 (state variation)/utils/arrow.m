%% global2local
function arrow(x, y, theta, color)

r = 0.5;                  % magnitude (length) of arrow to plot
u = r * cos(theta);
v = r * sin(theta);

q = quiver(x, y, u, v, color);
q.LineWidth = 1;
q.MaxHeadSize = 0.5;
q.AutoScale = 'off';

end