function [] = plot_xyz_time(x, y, z, t, namedPts)
figure;
szx = size(x);
n_plots = szx(1);
num_named_pts = size(namedPts);
names = string(["Satellite 1", "Satellite 2", "Satellite 3", "Satellite 4",...
    "Satellite 5", "Satellite 6", "Satellite 7", "Satellite 8"]);
satellites = string(["Satellite 1 ", "Satellite 2 ", "Satellite 3 ", "Satellite 4 ",...
    "Satellite 5", "Satellite 6", "Satellite 7", "Satellite 8"]);
directions = ['X', 'Y', 'Z'];
color = ['k', 'b', 'r', 'g', 'c', 'm', 'y'];
symbol = ['o', '+', 'd', 's', 'p'];

hold on
title('Spacecraft distance from satellite 1 vs. time')
xlabel('time (sec)')
ylabel('position (m)')
grid()

%plot3(0,0,0,strcat('b','x'), 'linewidth', 2,'DisplayName',names(1))
for i = 2:n_plots
%     plot(t, sqrt((x(i,:) - x(1,:)).^2 + (y(i,:) - y(1,:)).^2 + (y(i,:) - y(1,:)).^2),...
%         'linewidth', 2, 'DisplayName', satellites(i))
    plot(t, x(i,:) - x(1,:), 'linewidth', 2, 'DisplayName', strcat(satellites(i), directions(1)))
    plot(t, y(i,:) - y(1,:), 'linewidth', 2, 'DisplayName', strcat(satellites(i), directions(2)))
    plot(t, z(i,:) - z(1,:), 'linewidth', 2, 'DisplayName', strcat(satellites(i), directions(3)))
end

for j = 1:num_named_pts(2)
    plot(t(namedPts(j)), 0, strcat(color(1),symbol(j)),'linewidth', 2,...
        'DisplayName',strcat('Point', ' ',  int2str(j)));
end

legend()
hold off