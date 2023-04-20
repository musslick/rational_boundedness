clear all;
close all;
clc;
load('logfiles/Part2/GioManuscript_modifiedStroop_multiTraining_3P4F_4tasks_15_h100_ortho0');

y = mean(LCA_log);
y_sem = std(LCA_log)/sqrt(size(LCA_log,1));
x = 2:6;

color = [52 135 196]/255;

figure(1);
errorbar(x, y(x), y_sem(x), 'LineWidth', 3, 'Color', color);

xlabel(' ');
ylabel(' ');
ylim([0 1]);
xlim([1.5 6.5]);
set(gca, 'XTick', x);
set(gca, 'XTickLabel', 1:5);
set(gca, 'fontsize', 14);
