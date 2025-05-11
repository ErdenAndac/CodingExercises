clc;
clear;
close all;

% Load the data from deneme.m
rhoAir = 0.0023769;
rpm = 450;
omega = 450 * 2 * pi / 60;
tipR = 12;
rootR = 0.48;
numElem = 10;
dr = (tipR - rootR) / numElem;
discA = pi * tipR^2;
vi = sqrt(1390 / (2 * rhoAir * discA));
numBlade = 2;
theta = deg2rad(0);
beta = deg2rad(0);

% Calculate element positions
elemR = zeros(numElem,1);
for k = 1:numElem
    elemR(k) = rootR + (k-1) * dr + dr/2;
end

% Create body frame
xbody = zeros(numElem, 3);
ybody = zeros(numElem, 3);
zbody = zeros(numElem, 3);

for k = 1:numElem
    xbody(k,1) = 1;
    ybody(k,2) = 1;
    zbody(k,3) = 1;
end

% Initialize velocity components
Ubody = [0; 0; 0];  % Zero for hover condition

% Initialize matrix for storing local velocity components
% Dimensions: [numElem, numBlade, 3] where 3 represents [UR, UT, UP]
localVelocities = zeros(numElem, numBlade, 3);

% Time parameters
dt = 0.01;
azimuthInitial = zeros(1, numBlade);
angleAzimuth = 360/numBlade;
angle = 90;

for k = 1:numBlade
    azimuthInitial(1,k) = mod(deg2rad(angle + (k-1) * angleAzimuth), 2*pi);
end

% Plot setup
figure('Position', [100, 100, 1200, 800]);
hold on;
grid on;
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');

% Set axis limits
xlim([-12 12]);
ylim([-12 12]);
zlim([-12 12]);

% Enable zoom and rotate
zoom on;
rotate3d on;

% Plot body frame
quiver3(0,0,0,1,0,0,'r','LineWidth',2); % X-axis
quiver3(0,0,0,0,1,0,'g','LineWidth',2); % Y-axis
quiver3(0,0,0,0,0,1,'b','LineWidth',2); % Z-axis

% Number of time steps to show
numTimeSteps = 2;

% Plot blade elements for each time step
for i = 1:numTimeSteps
    time = dt * (i-1);
    
    % Clear previous plot except body frame
    cla;
    quiver3(0,0,0,1,0,0,'r','LineWidth',2);
    quiver3(0,0,0,0,1,0,'g','LineWidth',2);
    quiver3(0,0,0,0,0,1,'b','LineWidth',2);
    
    title(sprintf('Time Step: %d, Time: %.2f s', i, time));
    
    for j = 1:numBlade
        azimuth(j) = mod(azimuthInitial(j) + omega * time, 2*pi);
        % psi(j) = -((azimuthInitial(j) + omega * time) - deg2rad(angle));
        psi (j) = -azimuth(j);

        % Transformation matrices
        % z-axis (rotation of the blades around the hub)
        T1 = [cos(psi(j)),   sin(psi(j)),   0
              -sin(psi(j)),  cos(psi(j)),   0
              0,            0,              1];
        
        % y-axis (blade twist and pitch)
        T2 = [cos(theta),   0,  -sin(theta)
              0,            1,  0
              sin(theta),   0,   cos(theta)];

        % x-axis
        T3 = [1,    0,             0
              0,    cos(beta),     sin(beta)
              0,    -sin(beta),    cos(beta)];

        % Calculate hub frame (T2*T3 applied to body frame)
        xhub = zeros(numElem, 3);
        yhub = zeros(numElem, 3);
        zhub = zeros(numElem, 3);
        
        for k = 1:numElem
            xhub(k,:) = (T2 * T3) * xbody(k,:)';
            yhub(k,:) = (T2 * T3) * ybody(k,:)';
            zhub(k,:) = (T2 * T3) * -zbody(k,:)';
        end

        % Calculate blade frame (T1 applied to hub frame)
        xblade = zeros(numElem, 3);
        yblade = zeros(numElem, 3);
        zblade = zeros(numElem, 3);

        % Transform coordinates and calculate velocities
        for k = 1:numElem
            % Transform coordinates
            xblade(k,:) = T1 * xhub(k,:)';
            yblade(k,:) = T1 * yhub(k,:)';
            zblade(k,:) = T1 * zhub(k,:)';
            
            % Calculate position of each element
            pos = elemR(k) * [cos(azimuth(j)), sin(azimuth(j)), 0];
            
            % Calculate velocities
            Uhub = (T2 * T3) * Ubody;
            Ublade = T1 * Uhub;
            
            % Calculate local velocity components [UR UT UP]
            UR = Ublade(1);
            UT = -Ublade(2) + omega * elemR(k);
            UP = -Ublade(3) + vi;
            
            % Store local velocity components in the matrix
            localVelocities(k, j, 1) = UR;  % UR
            localVelocities(k, j, 2) = UT;  % UT
            localVelocities(k, j, 3) = UP;  % UP
            
            % Plot element position
            scatter3(pos(1), pos(2), pos(3), 'filled', 'MarkerFaceColor', [0.5 0.5 0.5]);
            
            % Plot hub frame at each element
            quiver3(pos(1), pos(2), pos(3), xhub(k,1), xhub(k,2), xhub(k,3), 'r--', 'LineWidth', 1);
            quiver3(pos(1), pos(2), pos(3), yhub(k,1), yhub(k,2), yhub(k,3), 'g--', 'LineWidth', 1);
            quiver3(pos(1), pos(2), pos(3), zhub(k,1), zhub(k,2), zhub(k,3), 'b--', 'LineWidth', 1);
            
            % Plot blade frame at each element
            quiver3(pos(1), pos(2), pos(3), xblade(k,1), xblade(k,2), xblade(k,3), 'r', 'LineWidth', 1);
            quiver3(pos(1), pos(2), pos(3), yblade(k,1), yblade(k,2), yblade(k,3), 'g', 'LineWidth', 1);
            quiver3(pos(1), pos(2), pos(3), zblade(k,1), zblade(k,2), zblade(k,3), 'b', 'LineWidth', 1);
        end
    end
    
    % Add legend
    legend('Body X', 'Body Y', 'Body Z', 'Blade Elements', 'Hub X', 'Hub Y', 'Hub Z', ...
        'Blade X', 'Blade Y', 'Blade Z');
    
    % Set view
    view(2);
    
    % Pause to see the animation
    pause(0.02);
end 

% Display the stored velocity components
disp('Local velocity components matrix:');
disp('Dimensions: [numElem, numBlade, 3] where 3 represents [UR, UT, UP]');
disp('First slice (UR components):');
disp(localVelocities(:,:,1));
disp('Second slice (UT components):');
disp(localVelocities(:,:,2));
disp('Third slice (UP components):');
disp(localVelocities(:,:,3));

rad2deg(azimuthInitial)
rad2deg(azimuth)
rad2deg(psi)