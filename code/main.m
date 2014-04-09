
%---------------------------------
% Initialization of the variables
function main
    global N N phi zeta theta_min theta_max delta alpha
    'start'
    N = 5; %Number of bees
    M = 2; %Number of tasks

    theta = ones(N,M)*500;
    x = ones(N,M)*0.1; 
    s  = ones(1,M);

    initial_condition = [reshape(theta, N*M,1); reshape(x, N*M,1);reshape(s,M,1)]; %intial condition for the ode_solver
    
    alpha = 3;
    delta = 1;
    p=0.2;
    zeta=10;
    phi=1;
    sigma=0.1;
    dt=1;
    
    theta_min=0; %we need to set theta_min
    theta_max=1000; %we need to set theta_max
    
    disp(dt)
    ode(1,initial_condition)
end


function dy = ode(t,y)
    global M N phi zeta theta_min theta_max delta alpha
    dy = zeros(2*N*M+M,1);
    
    %theta matrix
    for i=1:M
        for j=1:N
            factor = 0;
            y((i-1)*N + j, 1)
            if (y((i-1)*N + j, 1) > theta_min) && (y((i-1)*N + j, 1)<theta_max)
                factor = 1;
            end 
            
            dy((i-1)*N + j, 1) = [(1 - y((i-1)*N + j +N*M, 1))*phi- y((i-1)*N + j +N*M, 1)*zeta]*factor;
            %xij: dy((i-1)*N + j +N*M, 1)  
        end
    end
    
    %x matrix
    for i=1:M
        for j=1:N
            dy(M*N+(i-1)*N + j, 1) = 0.6;
        end
    end 
    
    %s matrix
    for i=1:M
        dy(2*M*N+i) = dy(2*M*N+i) + delta;
        for j=1:N
            dy(2*M*N+i)= dy(2*M*N+i) - alpha/N*y(M*N+(i-1)*N+j);
        end
    end

    dy
    
    
end



