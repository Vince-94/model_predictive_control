function [plotting] = plot_reference(x0, ref)

figure(1); hold on; grid on;
    title('Trajectory');
    plot(x0(1), x0(2), 'ok');
    plot(ref(1), ref(2), 'xk');


figure(2); hold on; grid on;
subplot(1,2,1); hold on; grid on;
    title('x');
%     plot(target.time, target.data(:,1), 'r');
subplot(1,2,2); hold on; grid on;
    title('y');
%     plot(target.time, target.data(:,2), 'r');
    
    
end