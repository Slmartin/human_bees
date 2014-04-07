
%---------------------------------
% Initialization of the variables
function main
    global N
    N = 5; %Number of bees
    global M 
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
    dt=1
    
    disp(dt)
    ode(1,1)
end


function dy = ode(t,y)
    global M N
    dy = zeros(2*N*M+M,1);
    %theta matrix
    for i=1:M
        for j=1:N
            dy((i-1)*N + j, 1) = 0.5
        end
    end
    %x matrix
    for i=1:M
        for j=1:N
            dy(M*N+(i-1)*N + j, 1) = 0.6
        end
    end 
    for i=1:M
        dy(2*M*N+i) = 0.7
    end
    M
    N
    
end


