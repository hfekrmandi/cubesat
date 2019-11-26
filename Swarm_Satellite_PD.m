clc
clear
close all
%% Input related things
delta_t = 0.01; %time stamp.
N = 3; % Number of Agents.
n_orbits = 1;
iteration = n_orbits * 530000;
%iteration = 100000;
%n_points = 100; % Very fast plot
n_points = 500; % Slow plot
% Orbittal period at 100km is 88.3 min; 53000 0.1 secs
format long
%% Constant Definitions
% PD controller gains
kp = -0.000001;
kd = -0.001;

mu = 3.986 * 10^14; % [m^3 / s^-2] gravitational parameter of primary body [earth in this case.]
earth_radius = 6.371 * 10^6; % [m] radius of primary body [earth in this case]
h = 100000; % [m] altitude of orbit from primary body
orbit_radius = h + earth_radius; % [m] radius of orbit from ceter of primary body
v = sqrt(mu / orbit_radius); % [m/s] Orbital Speed
omega = sqrt(mu / orbit_radius^3); % [rad/s] rotational speed.
T_earth = sqrt(4 * pi * orbit_radius^3 / mu); % [s] Orbital period

x = orbit_radius - sqrt(orbit_radius^2 + 1);

% Define position and velocity of satellites
spacecraft1_pos = [orbit_radius 0 0; 0 0 0];
spacecraft1_vel = [0 v 0; 0 0 0];
spacecraft2_pos_rel = [100 0 0; 0 0 0];
spacecraft2_vel_rel = [0 0 0; 0 0 0];
spacecraft3_pos_rel = [-100 0 100; 0 0 0];
spacecraft3_vel_rel = [0 0 0; 0 0 0];
spacecraft1_input = [0 0 0];
spacecraft2_input = [0 0 0];
spacecraft3_input = [0 0 0];

relation = [1,2; 1,3; 2,1; 2,3; 3,1; 3,2];

poses = [spacecraft1_pos;spacecraft1_pos + spacecraft2_pos_rel;...
    spacecraft1_pos + spacecraft3_pos_rel]; 
vels = [spacecraft1_vel;spacecraft1_vel + spacecraft2_vel_rel;...
    spacecraft1_vel + spacecraft3_vel_rel];
inputs = [spacecraft1_input;spacecraft2_input;spacecraft3_input];

%% Matrix and State Space Definitions
n_plots = 3;
x_list = zeros(n_plots, iteration);
y_list = zeros(n_plots, iteration);
z_list = zeros(n_plots, iteration);
u_list = zeros(18, iteration);
t_list = zeros(1, iteration);
x_max = 1;
y_max = 1;
z_max = 1;
x_min = 0;
y_min = 0;
z_min = 0;

zero = zeros(3);
I = eye(3);

% Hadaegh Matrix Definitions
D_0 = [3,0,0;0,0,0;0,0,-1];
S_0 = [0,2,0;-2,0,0;0,0,0];
A_0 = [zero,I; omega^2*D_0,omega*S_0];
A_expm = expm(A_0 * delta_t);
A_expm_kron = kron(eye(N), A_expm);

% Fixed Matrix Definitions
D_0 = [-3,0,0;0,-3,0;0,0,-1];
S_0 = [0,2,0;-2,0,0;0,0,0];
A_0 = [zero,I; omega^2*D_0,omega*S_0];
A_fixed = expm(A_0 * delta_t);
A_fixed_kron = kron(eye(N), A_fixed);
B_0 = [zero;I]; % inital B matrix
fun = @(tau)expm(A_0*(delta_t-tau))*B_0;
B_temp = integral(fun,0,delta_t,'ArrayValued',true);
B = kron(eye(N), B_0);
%A_fixed_kron = kron(eye(N), A_expm);

% Nonlinear Definitions
% A_Simple = zeros(6);
% A_Simple(1,4) = 1;
% A_Simple(2,5) = 1;
% A_Simple(3,6) = 1;
% A_Simple_kron = kron(eye(N), A_Simple);

% Turn initial pos and vel into state variable
X = [];
for i = 1:N
    X_pos = poses(2*i - 1,:);
    X_pos = reshape(X_pos, 1, []);
    X_vel = vels(2*i - 1,:);
    X_vel = reshape(X_vel, 1, []);
    X = [X; [X_pos, X_vel].'];
end
X_expm = X;
X_fixed = X;

for i=1:iteration
    % Compute nonlinear orbital dynamics
%     for j = 1:3
%         acc = -mu / (X(6*j - 5)^2 + X(6*j - 4)^2 + X(6*j - 3)^2);
%         acc_x = acc / sqrt(X(6*j - 5)^2 + X(6*j - 4)^2 + X(6*j - 3)^2);
%         acc_y = acc / sqrt(X(6*j - 5)^2 + X(6*j - 4)^2 + X(6*j - 3)^2);
%         acc_z = acc / sqrt(X(6*j - 5)^2 + X(6*j - 4)^2 + X(6*j - 3)^2);
%         A_Simple_kron(6*j - 2, 6*j - 5) = acc_x;
%         A_Simple_kron(6*j - 1, 6*j - 4) = acc_y;
%         A_Simple_kron(6*j - 0, 6*j - 3) = acc_y;
%     end
%     
%     X_dot = A_Simple_kron * X;
%     X_next = X + delta_t * X_dot;
%     X = X_next;

    % Conpute position error
    goal = [0 0 0 0 0 0 10 -10 0 0 0 0 -10 -10 0 0 0 0]';
    
    K = [zeros(18,18)];
    K(10,:) = [-kp 0 0 -kd 0 0 kp 0 0 kd 0 0 0 0 0 0 0 0];
    K(11,:) = [0 -kp 0 0 -kd 0 0 kp 0 0 kd 0 0 0 0 0 0 0];
    K(12,:) = [0 0 -kp 0 0 -kd 0 0 kp 0 0 kd 0 0 0 0 0 0];
    K(16,:) = [-kp 0 0 -kd 0 0 0 0 0 0 0 0 kp 0 0 kd 0 0];
    K(17,:) = [0 -kp 0 0 -kd 0 0 0 0 0 0 0 0 kp 0 0 kd 0];
    K(18,:) = [0 0 -kp 0 0 -kd 0 0 0 0 0 0 0 0 kp 0 0 kd];
    %K = [zeros(18,18)];
    
    % Compute Linear Fixed State
    X_next_fixed = A_fixed_kron * X_fixed + K * (X_fixed - goal);
    X_fixed = X_next_fixed;
    
    % Plot it
    for j = 1:3
        %x_list(j,i) = X(6*j - 5,1) - X(1,1);
        %y_list(j,i) = X(6*j - 4,1) - X(2,1);
        %z_list(j,i) = X(6*j - 3,1) - X(3,1);
        x_list(j,i) = X_fixed(6*j - 5,1);
        y_list(j,i) = X_fixed(6*j - 4,1);
        z_list(j,i) = X_fixed(6*j - 3,1);
        %u_list(:,i) = U;
        if (i > 1)
            t_list(i) = t_list(i - 1) + delta_t;
        end
    end
end

%% Plot relative positions
dispNames = string(["Satellite 1", "Satellite 2", "Satellite 3"]);
namedPts = [1 100000 200000 300000 400000];
plot_relative(x_list, y_list, z_list, namedPts, dispNames);
plot_absolute(x_list, y_list, z_list, namedPts, dispNames);
plot_xyz_time(x_list, y_list, z_list, t_list, namedPts);

%% plot control signals
% figure
% hold on
% grid()
% plot(t_list, u_list(10,:), 'DisplayName', 'Satellite 2 X')
% plot(t_list, u_list(11,:), 'DisplayName', 'Satellite 2 Y')
% plot(t_list, u_list(12,:), 'DisplayName', 'Satellite 2 Z')
% plot(t_list, u_list(16,:), 'DisplayName', 'Satellite 3 X')
% plot(t_list, u_list(17,:), 'DisplayName', 'Satellite 3 Y')
% plot(t_list, u_list(18,:), 'DisplayName', 'Satellite 3 Z')
% xlabel('time (seconds)')
% ylabel('Force (N)')
% legend()