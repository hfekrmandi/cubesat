% Plots the relative position of n satellites relative to the first

function [] = plot_relative(x, y, z, namedPts, names)
figure;
szx = size(x);
n_plots = szx(1);
num_named_pts = size(namedPts);
names = string(["Satellite 1", "Satellite 2", "Satellite 3", "Satellite 4",...
    "Satellite 5", "Satellite 6", "Satellite 7", "Satellite 8"]);
color = ['k', 'b', 'r', 'g', 'c', 'm', 'y'];
symbol = ['o', '+', 'd', 's', 'p'];

hold on
title('Spacecraft Positions Relative to Satellite 1')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid()

plot3(0,0,0,strcat('b','x'), 'linewidth', 2,'DisplayName',names(1))
for i = 2:n_plots
    plot3(x(i,:) - x(1,:),y(i,:) - y(1,:),z(i,:) - z(1,:),...
        'Color',color(mod(i - 1, length(color)) + 1),...
        'linestyle','-','LineWidth',1,'DisplayName',names(i))
end

for j = 1:num_named_pts(2)
    plot3([x(2,namedPts(j)) - x(1,namedPts(j)), x(3,namedPts(j)) - x(1,namedPts(j))],...
        [y(2,namedPts(j)) - y(1,namedPts(j)), y(3,namedPts(j)) - y(1,namedPts(j))],...
        [z(2,namedPts(j)) - z(1,namedPts(j)), z(3,namedPts(j)) - z(1,namedPts(j))],...
        strcat(color(1),symbol(j)),'linewidth', 2,...
        'DisplayName',strcat('Point', ' ',  int2str(j)));
end
legend()
hold off