function [] = plot_absolute(x, y, z, namedPts, names)
figure;
szx = size(x);
n_plots = szx(1);
num_named_pts = size(namedPts);
color = ['k', 'b', 'r', 'g', 'c', 'm', 'y'];
symbol = ['o', '+', 'd', 's', 'p'];

hold on
%plot3(0,0,0,'bo', 'linewidth', 2,'DisplayName','Earth Center')
for j = 1:num_named_pts(2)
    plot3(x(1,namedPts(j)),y(1,namedPts(j)),z(1,namedPts(j)),...
        strcat(color(1),symbol(j)),'linewidth', 2,...
        'DisplayName',strcat('Point', ' ', int2str(j)))
end

title('Absolute Spacecraft Positions')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
grid()

for i = 1:n_plots
    plot3(x(i,:),y(i,:),z(i,:), 'Color',color(mod(i - 1, length(color)) + 1),...
        'linestyle','-','LineWidth',1,'DisplayName',names(i))
end

[x,y,z] = sphere(100);
r = 6.371 * 10^6;
surf(x*r, y*r, z*r, 'FaceColor', [0 0 1], 'EdgeColor', 'none', 'DisplayName', 'Earth', 'FaceAlpha', 0.1);

legend()
hold off