% Robert Heeter
% BIOE 391 Numerical Methods
% HOMEWORK 10 MATLAB SCRIPT

clc, clf, clear, close all

%% P1. PROBLEM 23.27
disp('P1. PROBLEM 23.27');

tint = [0 20]; % time interval (days)
ivals = [0.05 5]; % inital values X(0) = 0.05, S(0) = 5

% PART A
tic;
[t,I] = ode23(@bacteria,tint,ivals); % use @bacteria function below for system of differential equations
time = toc; % execution timing
tXS = [t(:) I(:,1) I(:,2)]'; % combine results into one matrix

figure
hold on
plot(t,I(:,1),'-r','LineWidth',1.5);
plot(t,I(:,2),'-b','LineWidth',1.5);
hold off
grid on
xlabel('Time (t) [days]','FontSize',12,'FontWeight','bold');
ylabel('Bacterial biomass (X) or substrate concentration (S)','FontSize',12,'FontWeight','bold');
title('A. ode23 (rel tol = 1e-3)','FontSize',14,'FontWeight','bold');
legend('X (bacterial biomass)','S (substrate conc.)','FontSize',12,'FontWeight','bold');

fprintf('PART A: ode23: rel tol = 1e-3, incorrect (neg) results, tictoc = %fs\n  t         X        S\n',time);
fprintf('  %7.4f  %7.4f  %7.4f  \n',tXS);
fprintf('\n');

% PART B
opt = odeset('RelTol',1e-6); % set relative tolerance = 1e-6

tic;
[t,I] = ode23(@bacteria,tint,ivals,opt);
tXS = [t(:) I(:,1) I(:,2)]';
time = toc;

figure
hold on
plot(t,I(:,1),'-r','LineWidth',1.5);
plot(t,I(:,2),'-b','LineWidth',1.5);
hold off
grid on
xlabel('Time (t) [days]','FontSize',12,'FontWeight','bold');
ylabel('Bacterial biomass (X) or substrate concentration (S)','FontSize',12,'FontWeight','bold');
title('B. ode23 (rel tol = 1e-6)','FontSize',14,'FontWeight','bold');
legend('X (bacterial biomass)','S (substrate conc.)','FontSize',12,'FontWeight','bold');

fprintf('PART B: ode23: rel tol = 1e-6, correct (pos) results, tictoc = %fs\n  t         X        S\n',time);
fprintf('  %7.4f  %7.4f  %7.4f  \n',tXS(1:3,1:150));
fprintf('  ...       ...      ...\n\n');

% PART C
fprintf('PART C:\n');
tic;
[t,I] = ode45(@bacteria,tint,ivals,opt);
time1 = toc;
tXS1 = [t(:) I(:,1) I(:,2)]';
fprintf('ode45: rel tol = 1e-6, correct (pos) results, tictoc = %fs\n',time1);

figure
hold on
plot(t,I(:,1),'-r','LineWidth',1.5);
plot(t,I(:,2),'-b','LineWidth',1.5);
hold off
grid on
xlabel('Time (t) [days]','FontSize',12,'FontWeight','bold');
ylabel('Bacterial biomass (X) or substrate concentration (S)','FontSize',12,'FontWeight','bold');
title('C,i. ode45 (rel tol = 1e-6)','FontSize',14,'FontWeight','bold');
legend('X (bacterial biomass)','S (substrate conc.)','FontSize',12,'FontWeight','bold');

tic;
[t,I] = ode113(@bacteria,tint,ivals,opt);
time2 = toc;
tXS2 = [t(:) I(:,1) I(:,2)]';
fprintf('ode113: rel tol = 1e-6, correct (pos) results, tictoc = %fs\n',time2);

figure
hold on
plot(t,I(:,1),'-r','LineWidth',1.5);
plot(t,I(:,2),'-b','LineWidth',1.5);
hold off
grid on
xlabel('Time (t) [days]','FontSize',12,'FontWeight','bold');
ylabel('Bacterial biomass (X) or substrate concentration (S)','FontSize',12,'FontWeight','bold');
title('C,ii. ode113 (rel tol = 1e-6)','FontSize',14,'FontWeight','bold');
legend('X (bacterial biomass)','S (substrate conc.)','FontSize',12,'FontWeight','bold');

tic;
[t,I] = ode15s(@bacteria,tint,ivals,opt);
time3 = toc;
tXS3 = [t(:) I(:,1) I(:,2)]';
fprintf('ode15s: rel tol = 1e-6, correct (pos) results, tictoc = %fs\n',time3);

figure
hold on
plot(t,I(:,1),'-r','LineWidth',1.5);
plot(t,I(:,2),'-b','LineWidth',1.5);
hold off
grid on
xlabel('Time (t) [days]','FontSize',12,'FontWeight','bold');
ylabel('Bacterial biomass (X) or substrate concentration (S)','FontSize',12,'FontWeight','bold');
title('C,iii. ode15s (rel tol = 1e-6)','FontSize',14,'FontWeight','bold');
legend('X (bacterial biomass)','S (substrate conc.)','FontSize',12,'FontWeight','bold');

tic;
[t,I] = ode23s(@bacteria,tint,ivals,opt);
time4 = toc;
tXS4 = [t(:) I(:,1) I(:,2)]';
fprintf('ode23s: rel tol = 1e-6, correct (pos) results, tictoc = %fs\n',time4);

figure
hold on
plot(t,I(:,1),'-r','LineWidth',1.5);
plot(t,I(:,2),'-b','LineWidth',1.5);
hold off
grid on
xlabel('Time (t) [days]','FontSize',12,'FontWeight','bold');
ylabel('Bacterial biomass (X) or substrate concentration (S)','FontSize',12,'FontWeight','bold');
title('C,iv. ode23s (rel tol = 1e-6)','FontSize',14,'FontWeight','bold');
legend('X (bacterial biomass)','S (substrate conc.)','FontSize',12,'FontWeight','bold');

tic;
[t,I] = ode23t(@bacteria,tint,ivals,opt);
time5 = toc;
tXS5 = [t(:) I(:,1) I(:,2)]';
fprintf('ode23t: rel tol = 1e-6, correct (pos) results, tictoc = %fs\n',time5);

figure
hold on
plot(t,I(:,1),'-r','LineWidth',1.5);
plot(t,I(:,2),'-b','LineWidth',1.5);
hold off
grid on
xlabel('Time (t) [days]','FontSize',12,'FontWeight','bold');
ylabel('Bacterial biomass (X) or substrate concentration (S)','FontSize',12,'FontWeight','bold');
title('C,v. ode23t (rel tol = 1e-6)','FontSize',14,'FontWeight','bold');
legend('X (bacterial biomass)','S (substrate conc.)','FontSize',12,'FontWeight','bold');

tic;
[t,I] = ode23tb(@bacteria,tint,ivals,opt);
time6 = toc;
tXS6 = [t(:) I(:,1) I(:,2)]';
fprintf('ode23tb: rel tol = 1e-6, correct (pos) results, tictoc = %fs\n\n',time6);

figure
hold on
plot(t,I(:,1),'-r','LineWidth',1.5);
plot(t,I(:,2),'-b','LineWidth',1.5);
hold off
grid on
xlabel('Time (t) [days]','FontSize',12,'FontWeight','bold');
ylabel('Bacterial biomass (X) or substrate concentration (S)','FontSize',12,'FontWeight','bold');
title('C,vi. ode23tb (rel tol = 1e-6)','FontSize',14,'FontWeight','bold');
legend('X (bacterial biomass)','S (substrate conc.)','FontSize',12,'FontWeight','bold');


%% P2. PROBLEM 23.32
disp('P2. PROBLEM 23.32');

cdata = [85.3 66.6 60.6 56.1 49.1 45.3 41.9 37.8 33.7 34.4 35.1;
    16.9 18.7 24.1 20.9 18.9 19.9 20.6 13.9 19.1 14.5 15.4;
    4.7 7.9 20.1 22.8 32.5 37.7 42.4 47 50.5 52.3 51.3]'; % concentration data

kguess = [0.15 0.15 0.15 0.15]; % initial guess for constants k
k = fminsearch(@conc_SSR,kguess,[],cdata); % use fminsearch to optimize values of k with additional functions @conc_SSR and @concentration
fprintf('k constants:\nk_12 = %f\nk_21 = %f\nk_31 = %f\nk_32 = %f\n\n',k); % display results


%% P3. PROBLEM 24.11
disp('P3. PROBLEM 24.11');

D = 1.5e-6; % D coefficient (cm^2/s)
k = 5e-6; % k coefficient (s^-1)
L = 4; % length of tube (cm)

step = 0.2; % step size for finite difference
A_i = 0.1; % A(x=0) boundary concentration (M)
A_f = 0; % A(x=4) boundary concentration (M)

coefa = zeros(L/step - 1,L/step - 1); % preallocate matrix of coefficients
dim = size(coefa);
diag_1 = 1:(dim(1)+1):(dim(1)*dim(2)); % indices for middle diagonal
diag_2 = dim(1)+1:(dim(1)+1):(dim(1)*dim(2)); % indices for upper diagonal
diag_3 = 2:(dim(1)+1):(dim(1)*dim(2)); % indices for lower diagonal

coefa(diag_1) = 2 + ((k/D)*step^2); % set values for tridiagonal matrix
coefa(diag_2) = -1;
coefa(diag_3) = -1;

b = zeros(L/step - 1,1); % set values for constant column vector
b(1) = A_i;
b(end) = A_f;

A_int = gaussseidel(coefa,b); % use gaussseidel to solve coef*A_int = b for internal concentrations
A_comp = [A_i; A_int(:); A_f]; % complete concentration vector including boundary conditions
x = 0:step:L;

figure % plot results
hold on
fplot(@(x) (-4.5365247157e-8).*exp(1.82574185835.*x) + 0.100000045365.*exp(-1.82574185835.*x),[0,4],'-m','LineWidth',1.5); % analytical solution
plot(x,A_comp,'.k','MarkerSize',15);
hold off
grid on
xlabel('Distance along tube (x) [cm]','FontSize',12,'FontWeight','bold');
ylabel('Concentration of compound A [M]','FontSize',12,'FontWeight','bold');
title('Compound A diffusion along 4-cm tube','FontSize',14,'FontWeight','bold');
legend('Analytical solution','Centered FD solutions (step = 0.2cm)','FontSize',12,'FontWeight','bold');


%% P4. PROBLEM 24.13
disp('P4. PROBLEM 24.13');

D = 0.1; % constant (m^2/min)
U = 1; % constant (m/min)
k1 = 3; % constant (1/min)
k2 = 1; % constant (1/min)
L = 0.5; % length of reactor (m)

step = 0.05; % step size for finite difference
ca_i = 10; % input concentrations of A, B, & C
cb_i = 0;
cc_i = 0;

% PRODUCT A
coefa = zeros(L/step + 1,L/step + 1); % preallocate matrix of coefficients
dim = size(coefa);
diag_1 = 1:(dim(1)+1):(dim(1)*dim(2)); % indices for middle diagonal
diag_2 = dim(1)+1:(dim(1)+1):(dim(1)*dim(2)); % indices for upper diagonal
diag_3 = 2:(dim(1)+1):(dim(1)*dim(2)); % indices for lower diagonal

coefa(diag_1) = (-2*D/(step^2)) - k1; % set values for tridiagonal matrix
coefa(diag_2) = (D/(step^2)) - (U/(2*step));
coefa(diag_3) = (D/(step^2)) + (U/(2*step));

coefa(1,1:3) = [U+(3*D/(2*step)), -2*D/(step), D/(2*step)]; % boundary conditions
coefa(end,end-1:end) = [-2*D/(step^2), (2*D/(step^2))+k1];

b = zeros(L/step + 1,1); % set values for constant column vector
b(1) = U*ca_i; % boundary conditions
b(end) = 0;

ca_comp = gaussseidel(coefa,b); % use gaussseidel to solve coefa*ca_comp = b for all concentrations

% PRODUCT B
coefb = zeros(L/step + 1,L/step + 1); % preallocate matrix of coefficients
dim = size(coefb);
diag_1 = 1:(dim(1)+1):(dim(1)*dim(2)); % indices for middle diagonal
diag_2 = dim(1)+1:(dim(1)+1):(dim(1)*dim(2)); % indices for upper diagonal
diag_3 = 2:(dim(1)+1):(dim(1)*dim(2)); % indices for lower diagonal

coefb(diag_1) = (-2*D/(step^2)) - k2; % set values for tridiagonal matrix
coefb(diag_2) = (D/(step^2)) - (U/(2*step));
coefb(diag_3) = (D/(step^2)) + (U/(2*step));

coefb(1,1:3) = [U+(3*D/(2*step)), -2*D/(step), D/(2*step)]; % boundary conditions
coefb(end,end-1:end) = [-2*D/(step^2), (2*D/(step^2))+k2];

b = -1*k1*ca_comp; % set values for constant column vector
b(1) = U*cb_i; % boundary conditions
b(end) = k1*ca_comp(end);

cb_comp = gaussseidel(coefb,b); % use gaussseidel to solve coefb*cb_comp = b for all concentrations

% PRODUCT C
coefc = zeros(L/step + 1,L/step + 1); % preallocate matrix of coefficients
dim = size(coefc);
diag_1 = 1:(dim(1)+1):(dim(1)*dim(2)); % indices for middle diagonal
diag_2 = dim(1)+1:(dim(1)+1):(dim(1)*dim(2)); % indices for upper diagonal
diag_3 = 2:(dim(1)+1):(dim(1)*dim(2)); % indices for lower diagonal

coefc(diag_1) = (-2*D/(step^2)); % set values for tridiagonal matrix
coefc(diag_2) = (D/(step^2)) - (U/(2*step));
coefc(diag_3) = (D/(step^2)) + (U/(2*step));

coefc(1,1:3) = [U+(3*D/(2*step)), -2*D/(step), D/(2*step)]; % boundary conditions
coefc(end,end-1:end) = [-2*D/(step^2), (2*D/(step^2))];

b = -1*k2*cb_comp; % set values for constant column vector
b(1) = U*cc_i; % boundary conditions
b(end) = k2*cb_comp(end);

cc_comp = gaussseidel(coefc,b); % use gaussseidel to solve coefc*cc_comp = b for all concentrations
x = 0:step:L;

% PLOT RESULTS
figure
hold on
plot(x,ca_comp,'-r','LineWidth',1.5);
plot(x,cb_comp,'-g','LineWidth',1.5);
plot(x,cc_comp,'-b','LineWidth',1.5);
hold off
grid on
xlabel('Distance along reactor (x) [m]','FontSize',12,'FontWeight','bold');
ylabel('Concentration of compound [M]','FontSize',12,'FontWeight','bold');
title('Compound A, B, & C diffusion along 0.5m reactor','FontSize',14,'FontWeight','bold');
legend('Compound A','Compound B','Compound C','FontSize',12,'FontWeight','bold');


%% P5. PROBLEM 24.14
disp('P5. PROBLEM 24.14');

D = 0.8; % diffusion coefficient (cm^2/day)
Df = 0.64; % diffusion coefficient (cm^2/day)
k = 0.1; % first-order rate constant (1/day)
L = 0.008; % length of diffusion layer (cm)
Lf = 0.004; % length of biofilm (cm)

step = 0.001; % step size for finite difference (cm)
ca_i = 100; % bulk liquid concentrations for A (mol/L)

coefa = zeros((L+Lf)/step + 1,(L+Lf)/step + 1); % preallocate matrix of coefficients
dim = size(coefa);
diag_1 = 1:(dim(1)+1):(dim(1)*dim(2)); % indices for middle diagonal
diag_2 = dim(1)+1:(dim(1)+1):(dim(1)*dim(2)); % indices for upper diagonal
diag_3 = 2:(dim(1)+1):(dim(1)*dim(2)); % indices for lower diagonal

% 0 < X < L
coefa(diag_1(1:(L/step))) = -2*D/(step^2); % set values for tridiagonal matrix
coefa(diag_2(1:(L/step))) = D/(step^2);
coefa(diag_3(1:((L/step)-1))) = D/(step^2);

% L < X < (L + Lf)
coefa(diag_1(((L/step)+2):end)) = (-2*Df/(step^2)) - k;
coefa(diag_2(((L/step)+2):end)) = Df/(step^2);
coefa(diag_3(((L/step)+1):end)) = Df/(step^2);

% X = 0, L, L+Lf
coefa(1,1:2) = [(2*D/(step^2)) (-1*D/(step^2))]; % boundary/transition conditions
coefa(9,8:10) = [(D/step) ((-1*(D+Df)/step) - (k*step/2)) (Df/step)];
coefa(end,end-1:end) = [(-2*Df/(step^2)) ((2*Df/(step^2)) + k)];

b = zeros((L+Lf)/step + 1,1); % set values for constant column vector
b(1) = D*ca_i/(step^2); % boundary conditions
b(end) = 0;

ca_comp = gaussseidel(coefa,b); % use gaussseidel to solve coefa*ca_comp = b for all concentrations
x = 0:step:(L+Lf);

figure
plot(x,ca_comp,'-m','LineWidth',1.5);
grid on
xlabel('Distance (x) [cm]','FontSize',12,'FontWeight','bold');
ylabel('Concentration of compound A [M]','FontSize',12,'FontWeight','bold');
title('Compound A diffusion and reaction along diffusion layer and biofilm','FontSize',14,'FontWeight','bold');

xca_comp = [x(:) ca_comp(:)]'; % combine results into one matrix
fprintf(' x (cm)   c_a (M)\n');
fprintf('%7.4f  %8.4f\n',xca_comp);
fprintf('\n');


%% P6. PROBLEM 22.22
disp('P6. PROBLEM 24.22');

% ATTEMPTED BUT DID NOT FINISH
coefa = zeros(L/step - 1,L/step - 1); % preallocate matrix of coefficients
dim = size(coefa);
diag_1 = 1:(dim(1)+1):(dim(1)*dim(2)); % indices for middle diagonal
diag_2 = dim(1)+1:(dim(1)+1):(dim(1)*dim(2)); % indices for upper diagonal
diag_3 = 2:(dim(1)+1):(dim(1)*dim(2)); % indices for lower diagonal

coefa(diag_1) = 2 + ((k/D)*step^2); % set values for tridiagonal matrix
coefa(diag_2) = -1;
coefa(diag_3) = -1;

b = zeros(L/step - 1,1); % set values for constant column vector
b(1) = A_i;
b(end) = A_f;

A_int = gaussseidel(coefa,b); % use gaussseidel to solve coef*A_int = b for internal concentrations
A_comp = [A_i; A_int(:); A_f]; % complete concentration vector including boundary conditions
x = 0:step:L;


%% Additional Functions

function XS = bacteria(~,I)
% ABOUT: Differential equation system for P#1; X = I(1) = bacterial
% biomass, S = I(2) = substrate concentration

Y = 0.75; % yield coefficient
k_max = 0.3; % maximum baterial growth rate
k_s = 1e-4; % half saturation constant


XS = [Y*k_max*I(2)*I(1)/(k_s+I(2)); -1*k_max*I(2)*I(1)/(k_s+I(2))];

end


function Cs = concentrations(~,c,k)
% ABOUT: Differential equation system for P#2; c = vector of concentration
% of reactions 1-3, k = vector of constants for reactions 1-3

Cs = [-1*k(1)*c(1) + k(2)*c(2) + k(3)*c(3); k(1)*c(1) - k(2)*c(2) - k(4)*c(2); k(4)*c(2) - k(3)*c(3)];

end


function SSR = conc_SSR(k,cdata)
% ABOUT: Sums squares of discrepances between model predictions and data
% for P#2

tspan = [0 1 2 3 4 5 6 8 9 10 12 15]; % time interval (days)
ivals = [100 0 0]; % inital values c1(0) = 100, c2(0) = c3(0) = 0

[~,c] = ode45(@concentrations,tspan,ivals,[],k); % use ode45 to solve system of differential equations in function @concentration

R = (c(2:end,:)-cdata).^2; % calculate sum of squares of differences
SSR = sum(R,'all');

end


function x = gaussseidel(A,b,es,maxit)
% ABOUT: Gauss-Seidel method from textbook .m file.
% INPUTS: A = coefficient matrix; b = right side vector; es = relative
% error threshold (default = 0.00001%); maxit = max iterations (default =
% 50)
% OUTPUTS: x = solution vector

if nargin < 2
    error('At least 2 input arguments required.')
end
if nargin<4 || isempty(maxit)
    maxit=50;
end
if nargin<3 || isempty(es)
    es=0.00001;
end

[m,n] = size(A);
if m~=n
    error('Matrix A must be square.');
end

C = A;
x = zeros(1,n);
d = zeros(1,n);

for i = 1:n
  C(i,i) = 0;
end
x = x';
for i = 1:n
  C(i,1:n) = C(i,1:n)/A(i,i);
end
for i = 1:n
  d(i) = b(i)/A(i,i);
end

iter = 0;
ea = 100;
while (1)
  xold = x;
  for i = 1:n
    x(i) = d(i)-C(i,:)*x;
    if x(i) ~= 0
      ea = abs((x(i) - xold(i))/x(i)) * 100;
    end
  end
  iter = iter+1;
  if max(ea)<=es || iter >= maxit
      break
  end
end

end
