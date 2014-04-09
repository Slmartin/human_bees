
%---------------------------------
% Initialization of the variables
function main
    global N N phi zeta theta_min theta_max delta alpha p
    'start'
    N = 5; %Number of bees
    M = 2; %Number of tasks

    theta = ones(N,M)*500;
    x = ones(N,M)*0.1; 
    s  = ones(1,M);

    initial_condition = [reshape(theta, N*M,1); reshape(x, N*M,1);reshape(s,M,1)]; %intial condition for the ode_solver
    initial_condition
    
    alpha = 3;
    delta = 1;
    p=0.2;
    zeta=10;
    phi=1;
    sigma=0.1;
    dt=1;
    
    theta_min=1; 
    theta_max=1000; 
    
 
    ode(1,initial_condition)
end


function dy = ode(t,y)
    global M N phi zeta theta_min theta_max delta alpha p
    dy = zeros(2*N*M+M,1);
    
    %theta matrix
    for j=1:M
        for i=1:N
            factor = 0;
            y((j-1)*N + i, 1)
            if (y((j-1)*N + i, 1) > theta_min) && (y((j-1)*N + i, 1)<theta_max)
                factor = 1;
            end 
            
            dy((j-1)*N + i, 1) = [(1 - y((j-1)*N + i +N*M, 1))*phi- y((j-1)*N + i +N*M, 1)*zeta]*factor;
            %xij: dy((i-1)*N + j +N*M, 1)  
        end
    end
    
    %x matrix
    for j=1:M
        for i=1:N
            sum = 0;
            for w=1:M
                sum = sum+y(M*N+(w-1)*N + i);
            end
            gaussian_sp = 0;
            dy(M*N+(j-1)*N + i, 1) = [y(M*N*2+j)^2/(y(M*N*2+j)^2+y((j-1)*N + i)^2)]*[1-sum]-p*y(M*N+(j-1)*N + i, 1)+gaussian_sp;
        end
    end 
    
    %s matrix
    for j=1:M
        dy(2*M*N+j) = dy(2*M*N+j) + delta;
        for i=1:N
            dy(2*M*N+j)= dy(2*M*N+j) - alpha/N*y(M*N+(j-1)*N+i);
        end
    end

    dy
    
    
end



