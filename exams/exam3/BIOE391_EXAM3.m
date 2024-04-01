% Robert Heeter
% BIOE 391 Numerical Methods
% EXAM 3 MATLAB SCRIPT

clc, clf, clear, close all

%% PROBLEM 1, PART A
disp('PROBLEM 1A');

xint = [0 1.5]; % tube length (mm)
Ci = 30; % boundary conditions
Cf = 0;

z1 = 10; % guess z value for first two shots
z2 = 50;

[~,y1] = ode45(@linearrd,xint,[Ci z1]); % first shot
[~,y2] = ode45(@linearrd,xint,[Ci z2]); % second shot

za = z1 + (((z2-z1)/(y2(end,1)-y1(end,1)))*(Cf-y1(end,1))); % linear regression
[x3,y3] = ode45(@linearrd,xint,[Ci za]); % final shot

figure % plot results
plot(x3,y3(:,1),'-m','LineWidth',1.5)
xlabel('Tube position (x) [mm]','FontSize',12,'FontWeight','bold');
ylabel('Nutrient concentration (C) [mM]','FontSize',12,'FontWeight','bold');
title('PART A: Linear model, nutrient concentration vs. position','FontSize',14,'FontWeight','bold');
grid on


%% PROBLEM 1, PART B
disp('PROBLEM 1B');

xint = [0 1.5]; % tube length (mm)
Ci = 30; % boundary conditionsc

zguess = 0; % initial guess value for za
[za,r] = fzero(@res,zguess); % minimize z to fit boundary condition with fzero using res function below
[x,y] = ode45(@nonlinearrd,xint,[Ci za]); % final shot

figure % plot results
plot(x,y(:,1),'-m','LineWidth',1.5)
xlabel('Tube position (x) [mm]','FontSize',12,'FontWeight','bold');
ylabel('Nutrient concentration (C) [mM]','FontSize',12,'FontWeight','bold');
title('PART B: Nonlinear model, nutrient concentration vs. position','FontSize',14,'FontWeight','bold');
grid on


%% PROBLEM 2, PART A-D
disp('PROBLEM 2A-D');

% PART A-C: b = 6
b = 6; % dispersion term
dyediffusion(b); % use dyediffusion function below to solve and graph

% PART E: b = -4,-2,0,2,4
for b = -4:2:4
    dyediffusion(b);
end


%% PROBLEM 2, PART E
disp('PROBLEM 2E');

% No code necessary; see solution document


%% Additional Functions

function dy = linearrd(x,y)
% ABOUT: Linear system of ODEs to model the reaction and diffusion of
% nutrient in a chamber.

F = 0.8; % influx constant (mM/s)
x0 = 0.4; % influx constant (mm)
k = 3e-3; % consumption constant (1/s)
D = 2.5e-4; % diffusion constant (mm^2/s)

dy = [y(2);(1/D)*((k*y(1))-(F*(1/(x^2+x0^2))))]; % system

end


function dy = nonlinearrd(x,y)
% ABOUT: Nonlinear system of ODEs to model the reaction and diffusion of
% nutrient in a chamber.

F = 0.8; % influx constant (mM/s)
x0 = 0.4; % influx constant (mm)
k = 3e-3; % consumption constant (1/s)
D = 2.5e-4; % diffusion constant (mm^2/s)
KM = 2.5; % reaction constant (mM)

dy = [y(2);(1/D)*((k*y(1)/(1+(y(1)/KM)))-(F*(1/(x^2+x0^2))))]; % system

end


function r = res(za)
% ABOUT: Determines residual of nonlinear ODE at boundary for a given input
% za initial guess.

xint = [0 1.5]; % tube length (mm)
Ci = 30; % boundary conditions
Cf = 0;

[~,y] = ode45(@nonlinearrd,xint,[Ci,za]); % shot
r = y(end,1)-Cf; % residual

end


function [] = dyediffusion(b)
% ABOUT: Crank-Nicolson formulation to solve dye diffusion PDE with
% dispersion term b.

% Total number of nodes (chosen for good accuracy)
imax=201; % spatial index
nmax=200; % time index

% Initial step sizes
dx = 1/(imax-1);
dt = dx^2;

% Preallocate concentration matrix
C = zeros(imax,nmax+1);

% Initial conditions
C(:,1) = 0; % C(x,0) = 0
x = (0:dx:1)'; % space vector
t = 0; % time = 0

% Dirichlet boundary conditions
C(1,:) = 0; % C(0,t) = 0
C(end,:) = 1; % C(1,t) = 1

% Iterate for each time step
for n = 1:nmax
    
    % Dirichlet boundary conditions
    C(1,n+1) = 1;
    C(end,n+1) = 1;
    
    % Set values for tridiagonal matrix coef for coef*C(:,n+1)=d
    coef = zeros(imax,imax); % C(~,n+1)
    dim = size(coef);
    
    diag_mid = 1:(dim(1)+1):(dim(1)*dim(2)); % indices for middle diagonal
    diag_up = dim(1)+1:(dim(1)+1):(dim(1)*dim(2)); % indices for upper diagonal
    diag_low = 2:(dim(1)+1):(dim(1)*dim(2)); % indices for lower diagonal
    
    coef(diag_mid) = (-1/(dx^2)) + (-1/dt); % C(i,n+1)
    coef(diag_up) = (1/(2*dx^2)) + (b/(4*dx)); % C(i+1,n+1)
    coef(diag_low) = (1/(2*dx^2)) + (-b/(4*dx)); % C(i-1,n+1)
    
    coef(1,1:2) = [1,0]; % following Dirichlet boundary conditions
    coef(end,end-1:end) = [0,1];
    
    % Set values for d vector
    d = zeros(imax,1); % d(i-1/i/i+1,n)
    
    d1 = ((-1/(2*dx^2))+(b/(4*dx))) * C((1:end-2),n); % terms for d vector
    d2 = ((1/(dx^2))+(-1/dt)) * C((2:end-1),n);
    d3 = ((-1/(2*dx^2))+(-b/(4*dx))) * C((3:end),n);
    d(2:end-1) = d1+d2+d3;
    
    d(1) = 0; % following Dirichlet boundary conditions
    d(end) = 1;
    
    % Solve system using backslash and update concentration matrix
    C(:,n+1) = coef\d;
    
    % Update time and increase delta t by 10% each step
    t = t+dt;
    dt = 1.1*dt;
    
end

% Plot C vs. x for various values of t
figure
hold on
for i = 1:4:nmax/4
    plot(x,C(:,i),'Color',[1,0,0+(i/(nmax/4))],'LineWidth',1);
end
for i = (nmax/4)+1:4:nmax
    plot(x,C(:,i),'Color',[1-((i-(nmax/4))/(3*nmax/4)),0,1],'LineWidth',1);
end
title(['Dye diffusion in 2D chamber, b = ', num2str(b)],'FontSize',14','FontWeight','bold');
xlabel('Position in chamber (x)','FontSize',12,'FontWeight','bold');
ylabel('Concentration (C)','FontSize',12,'FontWeight','bold');
legend('initial state (t = 0)','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','steady state','FontSize',12,'FontWeight','bold','Location','northwest')
axis([-0.1,1.1,-0.1,1.1]);
hold off
grid on

end
