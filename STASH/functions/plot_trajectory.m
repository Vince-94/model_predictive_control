function [plotting] = plot_trajectory(x, color)

figure(1); hold on; grid on;
    title('Trajectory');
    plot(x.data(:,1), x.data(:,2), color);


figure(2); hold on; grid on;
subplot(1,2,1); hold on; grid on;
    title('x');
    plot(x.time, x.data(:,1), color);
subplot(1,2,2); hold on; grid on;
    title('y');
    plot(x.time, x.data(:,2), color);

    
end