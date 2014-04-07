%---------------------------------
% Initialization of the variables


N = 5; %Number of bees
M = 2; %Number of tasks

x = ones(N,M)*0.1;
theta = ones(N,M)*500;
s  = ones(1,M);

alpha = 3;
delta = 1;
p=0.2;
zeta=10;
phi=1;
sigma=0.1;
dt=1;
