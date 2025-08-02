%% 
% *MCP LOOP FOR LQR PROBLEM*
% 
% *OPEN LOOP*

import casadi.*
opti = casadi.Opti();

T = opti.variable(); %time horizon
N = 100; % no of sampling points
T = 10;
dt = T/N;

X = opti.variable(4, N+1) %states over time

U = opti.variable(2, N) % control inputs over time

x0 = [-5; -4; 0; 0] %initial state
xf = [pi/2; 0; 0; 0] %final state

% parameters
b1 = 180.0;
b2 = 45.0;
b3 = 23.5;
b4 = 25.0;
b5 = 122.5;
c1 = -25.0;
g1 = 784.8;
g2 = 245.3;

Q = 1e0*eye(4);
R = 1e0*eye(2);

%dynamics
f = @(X, U) [X(3);
    X(4);
    -(b5*U(1) - b3*U(2) - b5*g1*cos(X(1)) - b4*U(2)*cos(X(2)) + b3*g2*cos(X(1) + X(2)) - b5*g2*cos(X(1) + X(2)) + b4*g2*cos(X(1) + X(2))*cos(X(2)) + b3*c1*X(3)^2*sin(X(2)) + b5*c1*X(4)^2*sin(X(2)) + b4*c1*X(3)^2*cos(X(2))*sin(X(2)) + 2*b5*c1*X(3)*X(4)*sin(X(2)))/(b4^2*cos(X(2))^2 - b1*b5 + b3^2 - b2*b5*cos(X(2)) + 2*b3*b4*cos(X(2)));
    (b3*U(1) - b1*U(2) - b3*g1*cos(X(1)) - b2*U(2)*cos(X(2)) + b4*U(1)*cos(X(2)) + b1*g2*cos(X(1) + X(2)) - b3*g2*cos(X(1) + X(2)) + b2*g2*cos(X(1) + X(2))*cos(X(2)) - b4*g2*cos(X(1) + X(2))*cos(X(2)) - b4*g1*cos(X(1))*cos(X(2)) + b1*c1*X(3)^2*sin(X(2)) + b3*c1*X(4)^2*sin(X(2)) + b2*c1*X(3)^2*cos(X(2))*sin(X(2)) + b4*c1*X(4)^2*cos(X(2))*sin(X(2)) + 2*b3*c1*X(3)*X(4)*sin(X(2)) + 2*b4*c1*X(3)*X(4)*cos(X(2))*sin(X(2)))/(b4^2*cos(X(2))^2 - b1*b5 + b3^2 - b2*b5*cos(X(2)) + 2*b3*b4*cos(X(2)))]

% constraints
Umin = -1000;
Umax = 1000;
Qdotmin = -3*pi/2;
Qdotmax = 3*pi/2;
% opti.subject_to(Umin <= U <= Umax);
opti.subject_to(X(:,1)==x0);
opti.subject_to(X(:,N+1)==xf);
% opti.subject_to(Qdotmin <= X(3,:)<= Qdotmax)
% opti.subject_to(Qdotmin <= X(4,:)<= Qdotmax)


J=0;
for k =1:N
    %Runge-kutta integration
%     k1 = f(X(:, k), U(:, k)); 
%     k2 = f(X(:,k)+dt/2*k1, U(:,k));
%     k3 = f(X(:,k)+dt/2*k2, U(:,k));
%     k4 = f(X(:,k)+dt*k3, U(:,k));
%     x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4)
    % Euler integration
    x_next = X(:, k) + dt*f(X(:, k), U(:, k));
    % constraints
    opti.subject_to(X(:, k+1)==x_next);
    opti.subject_to(Umin <= U <= Umax);
    opti.subject_to(Qdotmin <= X(3,:)<= Qdotmax)
    opti.subject_to(Qdotmin <= X(4,:)<= Qdotmax)
    % cost function
    J = J + dt*(X(:, k+1)'*Q*X(:, k+1) + U(:, k)'*R*U(:, k));
end

opti.minimize(J);

opti.solver('ipopt')
sol = opti.solve;
%% 
% *FORWARD KINEMATICS*

l1 = 0.5;
l2 = 0.5;
q1_open = sol.value(X(1,:));
q2_open = sol.value(X(2,:));
v1_open = sol.value(X(3,:));
v2_open = sol.value(X(4,:));
u1_open = sol.value(U(1,:));
u2_open = sol.value(U(2,:));
% forward kinematics
x_open = l1*cos(q1_open) + l2*cos(q1_open + q2_open);
y_open = l1*sin(q1_open) + l2*sin(q1_open + q2_open);
%% 
% *MPC LOOP (CLOSED LOOP)*

import casadi.*
opti = casadi.Opti();

dt = 0.1;% sampling time

Nmpc = 100;% ocp iterations
T = Nmpc*dt; % time horizon
%constraints
Umin = -1000;
Umax = 1000;
Qdotmin = -3*pi/2;
Qdotmax = 3*pi/2;

x0 = [-5; -4; 0; 0] %initial state
xf = [pi/2; 0; 0; 0] %final state

N           =   10;
nx          =   4;
nu          =   2;
xsim           =   zeros(nx,N+1);
usim           =   zeros(nu,N);
xsim(:,1)      =   x0;
Xsol        =   repmat(x0,1,Nmpc+1);
Usol        =   repmat([0;0],1,Nmpc);
X_state = [];
U_state = [];
%MPC loop
for k = 1:10
    disp(k);
    x0 = xsim(:, k);
    [f, Xsol, Usol] = OCPloop(x0, xf, dt, opti, Nmpc, Xsol, Usol);
    usim(:, k) = Usol(:, 10);
    %xsim(:, k+1) = xsim(:, k) + dt*f(xsim(:,1), usim(:, 1));
    xsim(:, k+1) = Xsol(:, 10) + dt*f(Xsol(:, 10), Usol(:, 10));
    X_state = [X_state, Xsol(:, 1:10)];
    U_state = [U_state, Usol(:, 1:10)];
    Xsol = Xsol(:, 11:end);
    Usol = Usol(:, 11:end);
    Nmpc = Nmpc -10;%reducing the ocp horizon every iteration
end
%% 
% *FORWARD KINEMATICS*

X_state = [X_state, xsim(:, 11)]; % appending the last state
l1 = 0.5;
l2 = 0.5;
q1 = X_state(1,:);
q2 = X_state(2,:);
% forward kinematics
x = l1*cos(q1) + l2*cos(q1 + q2);
y = l1*sin(q1) + l2*sin(q1 + q2);

%% 
% *PLOTS*

figure(1);
clf(1);
plot(q1', q2');
grid on;
hold on;
plot(q1_open', q2_open');
xlabel("Q1");
ylabel("Q2");
legend('closed loop', 'open loop')
title('Q1-Q2 plane');
hold off;

figure(2);
clf(2);
plot(x', y');
grid on;
hold on;
plot(x_open', y_open')
xlabel("X");
ylabel("Y");
plot(x(1), y(1), 'r.', 'MarkerSize', 10);  % Start point
plot(x(end), y(end), 'g+', 'MarkerSize', 10);  % End point
legend('closed loop', 'open loop', 'start point', 'end point')
title('X-Y plane');
hold off; 

figure(3);
clf(3);
plot(xsim(3, :)', 'b');
hold on;
plot(xsim(4, :)', 'r');
hold on
plot(v1_open(:, 1:10:101)', 'y');
hold on
plot(v2_open(:, 1:10:101)', 'c');
grid on;
plot(1:10+1, Qdotmin*ones(10+1,1), '-.g');
plot(1:10+1, Qdotmax*ones(10+1,1), '-.g');
xlabel("Steps");
ylabel("Velocities");
legend('Q1dot', 'Q2dot', 'Q1dot open', 'Q2dot open');
xlim([1 12]);
title('Velocities');
hold off;

% Plot Inputs
figure(4);
clf(4);
plot(usim(1, :)', 'b');
hold on;
plot(usim(2, :)', 'r');
hold on,
plot(u1_open(:, 1:10:100)', 'y');
hold on,
plot(u2_open(:, 1:10:100)', 'c');
plot(1:11, Umin*ones(11,1), '-.g');
plot(1:11, Umax*ones(11,1), '-.g');
grid on;
xlabel("Steps");
ylabel("Inputs");
ylim([-1200 1200]);
xlim([0 11]);
legend('U1', 'U2');
title('Inputs');
%% 
% *OCP LOOP FUNCTION*

function [f, Xsol, Usol] = OCPloop(x0, xf, dt, opti, Nmpc, Xsol, Usol)
    
    b1 = 180.0;
    b2 = 45.0;
    b3 = 23.5;
    b4 = 25.0;
    b5 = 122.5;
    c1 = -25.0;
    g1 = 784.8;
    g2 = 245.3;
    
    Q = 1e0*eye(4);
    R = 1e0*eye(2);
        
    X = opti.variable(4, Nmpc+1); %states over time
    U = opti.variable(2, Nmpc); % control inputs over time

    f = @(X, U) [X(3);
    X(4);
    -(b5*U(1) - b3*U(2) - b5*g1*cos(X(1)) - b4*U(2)*cos(X(2)) + b3*g2*cos(X(1) + X(2)) - b5*g2*cos(X(1) + X(2)) + b4*g2*cos(X(1) + X(2))*cos(X(2)) + b3*c1*X(3)^2*sin(X(2)) + b5*c1*X(4)^2*sin(X(2)) + b4*c1*X(3)^2*cos(X(2))*sin(X(2)) + 2*b5*c1*X(3)*X(4)*sin(X(2)))/(b4^2*cos(X(2))^2 - b1*b5 + b3^2 - b2*b5*cos(X(2)) + 2*b3*b4*cos(X(2)));
    (b3*U(1) - b1*U(2) - b3*g1*cos(X(1)) - b2*U(2)*cos(X(2)) + b4*U(1)*cos(X(2)) + b1*g2*cos(X(1) + X(2)) - b3*g2*cos(X(1) + X(2)) + b2*g2*cos(X(1) + X(2))*cos(X(2)) - b4*g2*cos(X(1) + X(2))*cos(X(2)) - b4*g1*cos(X(1))*cos(X(2)) + b1*c1*X(3)^2*sin(X(2)) + b3*c1*X(4)^2*sin(X(2)) + b2*c1*X(3)^2*cos(X(2))*sin(X(2)) + b4*c1*X(4)^2*cos(X(2))*sin(X(2)) + 2*b3*c1*X(3)*X(4)*sin(X(2)) + 2*b4*c1*X(3)*X(4)*cos(X(2))*sin(X(2)))/(b4^2*cos(X(2))^2 - b1*b5 + b3^2 - b2*b5*cos(X(2)) + 2*b3*b4*cos(X(2)))];

    Umin = -1000;
    Umax = 1000;
    Qdotmin = -3*pi/2;
    Qdotmax = 3*pi/2;
    opti.subject_to(X(:,1)==x0);
    
    opti.set_initial(X,Xsol);
    opti.set_initial(U,Usol);

    J=0;
    for k =1:Nmpc
        %Runge-kutta integration
    %     k1 = f(X(:, k), U(:, k)); 
    %     k2 = f(X(:,k)+dt/2*k1, U(:,k));
    %     k3 = f(X(:,k)+dt/2*k2, U(:,k));
    %     k4 = f(X(:,k)+dt*k3, U(:,k));
    %     x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4)
        % Euler integration
        x_next = X(:, k) + dt*f(X(:, k), U(:, k));
        opti.subject_to(X(:, k+1)==x_next);
        opti.subject_to(Umin <= U <= Umax);
        opti.subject_to(Qdotmin <= X(3,:)<= Qdotmax)
        opti.subject_to(Qdotmin <= X(4,:)<= Qdotmax)
    
        J = J + dt*(X(:, k+1)'*Q*X(:, k+1) + U(:, k)'*R*U(:, k)); % terminal penalty
    end
    
    J = J + 10000*(X(:, end)-xf)'*(X(:, end)-xf);
    
    opti.minimize(J);
    
%     opti.set_initial(X,repmat(x0,1,Nmpc+1))
%     opti.set_initial(U,repmat([0;0],1,Nmpc));
    opti.solver('ipopt')
    sol = opti.solve;
    Xsol = sol.value(X);
    Usol = sol.value(U);
end