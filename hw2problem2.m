% RBE 501 - Robot Dynamics - Fall 2021
% Homework 2, Problem 2
% Worcester Polytechnic Institute
%
% Instructor: L. Fichera <lfichera@wpi.edu>
% Last modified: 09/22/2021
clear, clc, close all
addpath('utils');

plotOn = true;
nTests = 20; % number of random test configurations

%% Create the manipulator
% robot length values (meters)
L0 = 0.3;
L1 = 0.3;
L2 = 0.3;

robot = SerialLink([Revolute('a', 0, 'd', L0, 'alpha', -pi/2, 'offset', pi/2), ...
                    Revolute('a', L1, 'd', 0, 'alpha', 0), ...
                    Revolute('a', L2, 'd', 0, 'alpha', pi/2), ...
                    Revolute('a', 0, 'd', 0, 'alpha', -pi/2), ...
                    Revolute('a', 0, 'd', 0, 'alpha', -pi/2, 'offset', -pi/2), ...
                    Revolute('a', 0, 'd', 0, 'alpha', 0)], 'name', 'Elbow Manipulator'); 

% Joint limits
qlim = [-pi/2  pi/2;  % q(1)
        -pi/4  pi/2;  % q(2)
        0      pi/3;  % q(3)
        -pi/2  pi/2;  % q(4)
        -pi/2  pi/2;  % q(5)
        -pi/2  pi/2]; % q(6)

% Display the manipulator in the home configuration
q = zeros(1,6);
robot.teach(q);

%% Part A - Calculating the screw axes
% Screw axis of all the joints
S = [0 0 1 0 0 0; -1 0 0 0 -L0 0; -1 0 0 0 -L0 L1; 0 0 1 (L1+L2) 0 0; -1 0 0 0 -L0 (L1+L2); 0 1 0 -L0 0 0]';


%% Part B - Calculating the forward kinematics with the Product of Exponentials formula
M = [0 1 0 0; 0 0 1 (L1+L2); 1 0 0 L0; 0 0 0 1]; % home configuration

fprintf('---------------------Forward Kinematics Test---------------------\n');
fprintf(['Testing ' num2str(nTests) ' random configurations.\n']);
fprintf('Progress: ');
nbytes = fprintf('0%%'); 
 
% Testing the forward kinematics for 20 random sets of joint variables
for ii = 1 : nTests
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('%0.f%%', ceil(ii/nTests*100));
    
    %Generating a random configuration
    q = [qlim(1,1) + (qlim(1,2) - qlim(1,1)) * rand(), ...
         qlim(2,1) + (qlim(2,2) - qlim(2,1)) * rand(), ...
         qlim(3,1) + (qlim(3,2) - qlim(3,1)) * rand(), ...
         qlim(4,1) + (qlim(4,2) - qlim(4,1)) * rand(), ...
         qlim(5,1) + (qlim(5,2) - qlim(5,1)) * rand(), ...
         qlim(6,1) + (qlim(6,2) - qlim(6,1)) * rand()];
    
    %Calculating the forward kinematics
    T = fkine(S,M,q);
    
    if plotOn
        robot.teach(q);
        title('Forward Kinematics Test');
    end
    
    %Checking the correctnedd of the fkine function developed by me
    assert(all(all(abs(double(robot.fkine(q)) - T) < 1e-10)));
end
 
fprintf('\nTest passed successfully.\n');

%% Part C - Calculating the Space Jacobian of the manipulator

fprintf('-------------------Differential Kinematics Test------------------\n');
fprintf(['Testing ' num2str(nTests) ' random configurations.\n']);
fprintf('Progress: ');
nbytes = fprintf('0%%'); 

% Testing the correctness of the Jacobian for 100 random sets of joint variables
for ii = 1 : nTests
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('%0.f%%', ceil(ii/nTests*100));
    
    % Generating a random configuration
    q = [qlim(1,1) + (qlim(1,2) - qlim(1,1)) * rand(), ...
         qlim(2,1) + (qlim(2,2) - qlim(2,1)) * rand(), ...
         qlim(3,1) + (qlim(3,2) - qlim(3,1)) * rand(), ...
         qlim(4,1) + (qlim(4,2) - qlim(4,1)) * rand(), ...
         qlim(5,1) + (qlim(5,2) - qlim(5,1)) * rand(), ...
         qlim(6,1) + (qlim(6,2) - qlim(6,1)) * rand()];
    
    % Calculating the Forward Kinematics
    T = fkine(S,M,q);
    
    % Calculating the Jacobian
    J = jacob0(S,q);
    
    if plotOn
        robot.teach(q);
        title('Differential Kinematics Test');
    end
    
    % Testing the correctness of the Jacob0 function dveeloped by me
    Jcoords = [-skew(T(1:3,4))*J(1:3,:)+J(4:6,:); J(1:3,:)];
    assert(all(all(abs(double(robot.jacob0(q)) - Jcoords) < 1e-10)));
end

fprintf('\nTest passed successfully.\n');

%% Part D - Inverse Kinematics
fprintf('----------------------Inverse Kinematics Test--------------------\n');
fprintf(['Testing ' num2str(nTests) ' random configurations.\n']);
fprintf('Progress: ');
nbytes = fprintf('0%%');

% Calculating the twist representing the robot's home pose
currentPose = MatrixLog6(M);
currentPose = [currentPose(3,2) currentPose(1,3) currentPose(2,1) currentPose(1:3,4)']';

% Initializing the current joint variables
currentQ = zeros(1,6);

if plotOn
    robot.teach(currentQ);
    h = triad('matrix', M, 'tag', 'Target Pose', 'linewidth', 2.5, 'scale', 0.5);
end
     
% Generating the test configurations
q = [linspace(0,pi/5,nTests);
     linspace(0,pi/8,nTests);
     linspace(0,pi/8,nTests);
     linspace(0,pi/4,nTests);
     linspace(0,pi/8,nTests);
     linspace(0,pi/6,nTests)];
 
for ii = 1 : nTests
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('%0.f%%', ceil(ii/nTests*100));
    
    % Generating the robot's target pose
    T = fkine(S,M,q(:,ii)');
    targetPose = MatrixLog6(T);
    targetPose = [targetPose(3,2) targetPose(1,3) targetPose(2,1) targetPose(1:3,4)']';
    
    if plotOn
        set(h, 'matrix', T);
        title('Inverse Kinematics Test');
        drawnow;
    end
    
    iterations = 0;
    x_graph = [];
    y_graph = [];
    
    % Inverse Kinematics
    while (norm(targetPose - currentPose) > 1e-3)
        J = jacob0(S,currentQ);
        lambda = 2.0;
        J_star = J'*pinv(J*J' + (lambda^2)*eye(6));
        deltaQ = J_star*(targetPose - currentPose);
        %deltaQ = pinv(J)*(targetPose - currentPose);
        
        iterations = iterations + 1;
        x_graph(iterations) = iterations; 
        y_graph(iterations) = (norm(targetPose - currentPose));
        error = (norm(targetPose - currentPose))
        
        figure(2)
        plot(x_graph,y_graph)
        title(sprintf('Damped-Least-Square-NEWTON-RAPHSON of Configuration: %d', ii))
        xlabel('iterations')
        ylabel('norm(target-pose - current-pose)')
        hold on
        
        currentQ = currentQ + deltaQ';
        
        T = fkine(S,M,currentQ);
        currentPose = MatrixLog6(T);
        currentPose = [currentPose(3,2) ...
                       currentPose(1,3) ...
                       currentPose(2,1) ...
                       currentPose(1:3,4)']';
        if plotOn
            try
                figure(1)
                robot.teach(currentQ);
                drawnow;
            catch e
                continue;
            end
        end
    end
    clf(figure(2));
end
fprintf('\nTest passed successfully.\n');
