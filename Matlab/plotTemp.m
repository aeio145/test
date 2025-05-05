%Plottemp vs position

figure();

plot(Position,Temp,'LineWidth',2);
hold on
grid on
xlabel("Position x (m)");
ylabel("Temperature (ºC)");
title("Temperature vs. Position");
xlim([min(Position) max(Position)]);
ylim([min(Temp) max(Temp)]);
set(gca,'XMinorTick', 'on', 'YMinorTick','on', 'XMinorGrid', 'on', 'YMinorGrid','on')

if save
    % Save the plot:
    set(gcf, 'Position', [0, 0, 1280, 720]); % Set the figure size.
    saveas(gcf, fullfile('Results', 'plotTemp.fig'));  % Save the plot as a FIG file.
    saveas(gcf, fullfile('Results', 'plotTemp.png'));  % Save the plot as a PNG file.
end

disp(">> Temperature plot completed.");
fprintf('\n');

