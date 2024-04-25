close all; clear; clc;
% Constants
T = 293.15; %K
k_B = 1.38*10^-23; %J/K
mu_0 = 4*pi*10^-7; %Tm/A

% Ink parameters
a = 7*10^-9; %7nm
chi = 16;
eta = 0.08; %Pa*s
gamma = 6*pi*eta*a;
N = 150;

% Wall definitions taken from printing bed
x_wall = [-0.5 0.5]*10.^-7;
y_wall = [-0.5 0.5]*10.^-7;
z_wall = [-0.5 0.5]*10.^-7;

% Magnetic Field
B = [400*10^-3 0 0 ]; %20 mT
H_0 = B/(mu_0*(1+chi));
m = 4*pi*(a.^3)*chi.*H_0/3;

% Magnetic Coupling Parameter
lambda = mu_0*pi*(a^3)*(chi^2)*(norm(H_0)^2)/(9*k_B*T);
phi_3d = (4*pi*a^3/3)*N/10^-21;


% Normalisation
F0 = (3*mu_0*norm(m)^2)/(4*pi*(2*a)^4);
tc = (12*pi*eta*a^2)/F0;
tau = 0.05;
dt = tau*tc;

% Initial Positions: A voxel from -1mu_m to 1mu_m in every direction
% Scale the starting positions by non-dimesionalising
vox_start = (-0.5*10^-7);
vox_end = (0.5*10^-7);

pos = zeros(N,3);
% Generate random points
i = 1;
while i <= N
    rnd = (vox_start + (vox_end-vox_start).*rand(1,3,'double'));
    pos(i,:) = rnd;
    for j = i-1:-1:1
        check = norm(rnd - pos(j,:));
        if check <(2*a)
            % Regenerate i
            i = i-1;
            break;
        end
    end
    i =i+1;
end



% Plotting
% First define a sphere
figure
set(gcf, 'Position',  [100, 100, 1100, 1200]);
[X_sp, Y_sp, Z_sp] = sphere;
X_sp = X_sp*a;
Y_sp = Y_sp*a;
Z_sp = Z_sp*a;
% Surface plot each sphere individually
sp = [];
colormap("sky")
%figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:N
    plot_sphere = surf(X_sp+pos(i,1), Y_sp+pos(i,2), Z_sp+pos(i,3));
    hold on;
    sp = [sp, plot_sphere];
end
xlim([-1*10^-7 1*10^-7])
ylim([-1*10^-7 1*10^-7])
zlim([-1*10^-7 1*10^-7])
xlabel('x');
ylabel('y');
zlabel('z');
annstr = sprintf('Magnetic Field = <%0.3f, %0.3f, %0.3f> T', B(1), B(2), B(3) ); % annotation text
annpos = [0.39 0.87 0.1 0.1]; % annotation position in figure coordinates
ha = annotation('textbox',annpos,'string',annstr);
ha.FontSize = 20;
annstr = sprintf('time = %f s',0); % annotation text
annpos = [0.8 0 0.1 0.1]; % annotation position in figure coordinates
ha = annotation('textbox',annpos,'string',annstr);
ha.FontSize = 15;
drawnow;

hold on;


for time = 0:dt:10
    Force = zeros(N,3);
    for i = 1:N-1
        for j = (i+1):N
            % In this loop we calculate forces between particles
            % Calculate distance between the two particles
            r = pos(j,:) - pos(i,:);

            % theta is the angle between m and r
            cos_theta = m*r'/(norm(m)*norm(r));
            theta = acos(cos_theta);

            % Unit vectors
            r_hat = r./norm(r);
            m_hat = m./norm(m);
            % Need to find unit vector in theta direction

            % Magnetic force
            mx = m_hat(1);
            my = m_hat(2);
            mz = m_hat(3);
            rx = r_hat(1);
            ry = r_hat(2);
            rz = r_hat(3);

            tx=(rx*sin(theta)*my^2 - mx*ry*sin(theta)*my + rx*sin(theta)*mz^2 - mx*rz*sin(theta)*mz)/(mx^2*ry^2 + mx^2*rz^2 - 2*mx*my*rx*ry - 2*mx*mz*rx*rz + my^2*rx^2 + my^2*rz^2 - 2*my*mz*ry*rz + mz^2*rx^2 + mz^2*ry^2);
            ty=(sin(theta)*(ry*mx^2 - my*rx*mx + ry*mz^2 - my*rz*mz))/(mx^2*ry^2 + mx^2*rz^2 - 2*mx*my*rx*ry - 2*mx*mz*rx*rz + my^2*rx^2 + my^2*rz^2 - 2*my*mz*ry*rz + mz^2*rx^2 + mz^2*ry^2);
            tz=(sin(theta)*(rz*mx^2 - mz*rx*mx + rz*my^2 - mz*ry*my))/(mx^2*ry^2 + mx^2*rz^2 - 2*mx*my*rx*ry - 2*mx*mz*rx*rz + my^2*rx^2 + my^2*rz^2 - 2*my*mz*ry*rz + mz^2*rx^2 + mz^2*ry^2);
            t_hat = [tx, ty, tz];

            F_r = (F0*(2*a/norm(r))^4)*(3*(cos_theta^2) - 1);
            F_theta = (F0*(2*a/norm(r))^4)*sin(2*theta);

            % Short Range Repulsive force between particles
            F_pp = (2*F0)*exp(-10*(norm(r)/(2*a) -1));

            if(norm(r)<(2*a))
                F_theta = 0;
                F_r = F0;
                F_pp = (2*F0);
            end

            F = (F_r.*r_hat)+(F_pp.*-1.*r_hat)+(F_theta.*t_hat);

            % Sum of particle-particle forces
            Force(i,:) = Force(i,:) + F;
            Force(j,:) = Force(j,:) - F;
        end

        % Now the forces that are not between particles but on the
        % particle individually

        % Short Range Repulsive force from wall to particle i
        dist_wall_x = pos(i,1) - x_wall;
        F_pw_x = sum((2*F0)*exp(-10.*abs(dist_wall_x)./(2.*a) -0.5));
        if (dist_wall_x(1)<0 || dist_wall_x(1)<(a)) %too close/crossed
            F_pw_x = (2*F0)*tau;
        elseif (dist_wall_x(2)>0 || abs(dist_wall_x(2))<(a)) 
            F_pw_x = -(2*F0)*tau;
        end

        dist_wall_y = pos(i,2) - y_wall;
        F_pw_y = sum((2*F0)*exp(-10.*abs(dist_wall_y)./(2.*a) -0.5));
        if (dist_wall_y(1)<0 || dist_wall_y(1)<(a)) %too close/crossed
            F_pw_y = (2*F0)*tau;
        elseif (dist_wall_y(2)>0 || abs(dist_wall_y(2))<(a)) 
            F_pw_y = -(2*F0)*tau;
        end


        dist_wall_z = pos(i,3) - z_wall;
        F_pw_z = sum((2*F0)*exp(-10.*abs(dist_wall_z)./(2.*a) -0.5));
        if (dist_wall_z(1)<0 || dist_wall_z(1)<(a)) %too close/crossed
            F_pw_z = (2*F0)*tau;
        elseif (dist_wall_z(2)>0 || abs(dist_wall_z(2))<(a)) 
            F_pw_z = -(2*F0)*tau;
        end


        % Add random brownian motion to particle i
        R = normrnd(0,1,[1,3])*12*pi*a*eta*k_B*T/dt;

        % Sum of Forces on i
        Force(i,:) = Force(i,:) + R + [F_pw_x F_pw_y F_pw_z];
    end


    % ALl forces on all particles have been calculated now time for solving
    % the differential equation
    % Everything is currently dimensional
    X = pos./(2*a); % non-dimensional position
    f = Force./F0; % non-dimensional force
    % Solve linear differntial equation
    pos_nd = X + f.*tau;
    % Convert to dimensional
    pos = pos_nd.*2.*a;

    %Plot
    delete(sp);
    delete(ha)
    sp = [];
    for i = 1:N
        plot_sphere = surf(X_sp+pos(i,1), Y_sp+pos(i,2), Z_sp+pos(i,3));
        hold on;
        sp = [sp, plot_sphere];
    end
    hold on;
    annstr = sprintf('time = %f s',time); % annotation text
    annpos = [0.8 0 0.1 0.1]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr);
    ha.FontSize = 15;
    drawnow;
end