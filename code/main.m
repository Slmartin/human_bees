

%---------------------------------
% Initialization of the variables
function main
    global M N phi zeta theta_min theta_max delta alpha p sigma
    'start'
    N = 10; %Number of bees
    M = 5; %Number of tasks

    theta = ones(N,M)*500;
    x = ones(N,M)*0.1;
    s = ones(1,M);

    initial_condition = [reshape(theta, N*M,1); reshape(x, N*M,1);reshape(s,M,1)]; %intial condition for the ode_solver
    
    alpha = M+1;  %This factor has been changed in order to get convergence
                  %We could also change the delta!
                  %This change has been made because it doesn't make much sense
                  %to get an infinitely increasing stimuli for all tasks.

    delta = 1; 
    p=0.2;
    zeta=10;
    phi=1;
    sigma=0.1;
    dt=1;
    end_time = 5000;
    
    theta_min=1;
    theta_max=1000;
    
    %[T,Y]=ode45(@ode_function, [0 3000], initial_condition);
    [T,Y]=euler_method(@ode_function, [0 end_time], initial_condition,0.5);
    

    %plot(T,Y(:,1:M*N),'-o')
    %plot(T,Y(:,13),'-r')
    subplot(2,2,1)
    plot(T,Y(1:M*N,:),'-')
    axis([0 end_time 0 1100])
    title('theta')
    
    subplot(2,2,2)
    plot(T,Y(M*N+1:M*N*2,:),'-')
    axis([0 end_time 0 1.5])
    title('x')
    
    subplot(2,2,3)
    plot(T,Y(M*N*2+1:M*N*2+M,:),'-')
    title('s')
    axis([0 end_time 0 400])
    
    %here is a first version of welfare implementation
    %of course we can use another function to estimate welfare

    welfare = 1 ./ (Y(M*N*2+1,:) + 1);
    for i=2:M
        welfare = welfare + 1./ (Y(M*N*2+i,:) + 1);
    end
    
    subplot(2,2,4)
    plot(T,welfare,'-')
    axis([0 end_time 0 3])
    title('welfare')
    

end


function dy = ode_function(t,y)
    global M N phi zeta theta_min theta_max delta alpha p sigma
    dy = zeros(2*N*M+M,1);
    
    %theta matrix
    for j=1:M
        for i=1:N
            factor = 0;
            y((j-1)*N + i, 1);
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
            %gaussian_sp = 0;
            %here is a attempt to get the gaussian centerd process
            %unfortunately the gaussian random number generator seems to be
            %way to slow!!!
            %%%gaussian_sp = normrnd(0,sigma);
            gaussian_sp = sigma*randn;
            
            dy(M*N+(j-1)*N + i, 1) = [y(M*N*2+j)^2/(y(M*N*2+j)^2+y((j-1)*N + i)^2)]*[1-sum]-p*y(M*N+(j-1)*N + i, 1)+gaussian_sp;
        end
    end 
  
    %s matrix
    for j=1:M
        %it makes no sense to have a negative value for the stiumuli
        sum=0;
        for i=1:N
            %it makes no sense to have a negative value for the stiumuli
            sum = sum + alpha/N*y(M*N+(j-1)*N+i);
        end
        
        if ~((y(2*M*N+j)<=0) && ((delta - sum)<0))
            dy(2*M*N+j)= dy(2*M*N+j) + delta - sum;
        end
    end
    
    
end


function [T,Y] = euler_method(fun,time, y0, h)

    T = [];
    Y = [];
    T = [T,time(1)];
    Y = [Y,y0];
    
    
    y_eul = y0;
    t_eul = time(1);
    
    number = (time(2)-time(1))/h;
    
    for j=1:number
        y_eul = y_eul + h*fun(t_eul,y_eul);
        t_eul = t_eul + h ;
        
        Y = [Y, y_eul];
        T = [T, t_eul];
    end
    

end 
