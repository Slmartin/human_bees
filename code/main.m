%---------------------------------
% Initialization of the variables

%Set seed for random number generation
RandStream.setGlobalStream(RandStream('mt19937ar','seed',42))
RandStream.getGlobalStream()

function main
    global M N phi zeta theta_min theta_max delta alpha p sigma beta
    'start'
    N = 5;                     %Number of bees
    M = 2;                      %Number of tasks

    theta = ones(N,M)*500;      %response treshold
    x = ones(N,M)*0.1;          %fraction of time spent by individual i on task j
    s = ones(1,M);              %stimulus intensity for task j

    initial_condition = [reshape(theta, N*M,1); reshape(x, N*M,1);reshape(s,M,1)]; %intial condition for the ode_solver
    
    alpha = M+1;                %This factor has been changed in order to get convergence
                                %We could also change the delta!
                                %This change has been made because it doesn't make much sense
                                %to get an infinitely increasing stimuli for all tasks.

    delta = 1;                  %describes increase of stimuli
    p=0.2;                      %Probability individual i gives up performing task j in time interval
    zeta=5;                     %learning %4
    phi=1;                      %forgetting
    sigma=0.07;                  %sigma for gaussian distribution
    beta=50;                    %prefactor for growth rate

    dt=0.9;
    end_time = 6000;

    %dt=0.01;
    %end_time = 5000;

    
    theta_min=1;
    theta_max=1000;
    
    %[T,Y]=ode45(@ode_function, [0 3000], initial_condition);
    [T,Y]=euler_method(@ode_function, [0 end_time], initial_condition,dt);
    

    %plot(T,Y(:,1:M*N),'-o')
    %plot(T,Y(:,13),'-r')
    %subplot(3,2,1)
    %plot(T,Y(1:M*N,:),'-')
    %axis([0 end_time 0 1100])
    %title('theta')
    
    %subplot(3,2,2)
    %plot(T,Y(M*N+1:M*N*2,:),'-')
    %axis([0 end_time -0.1 1.1])
    %title('x')
    
    %subplot(3,2,3)
    %plot(T,Y(M*N*2+1:M*N*2+M,:),'-')
    %title('s')
    %axis([0 end_time 0 400])
    
    %%%%%%%%%%%%%%%% Plotting section report %%%%%%%%%%%

 set(0,'DefaultAxesColorOrder',jet(10))
     subplot(1,3,1)
    set(gca,'fontsize',28)
    plot(T,Y(1:M*N,:),'-','LineWidth',1.5)
    axis([0 end_time 0 1050])
    xlabel('Time','FontSize',30)
    ylabel('\theta','FontSize',30)
     hleg1=legend('\theta_{11}','\theta_{21}','\theta_{31}','\theta_{41}','\theta_{51}','\theta_{12}','\theta_{22}','\theta_{32}','\theta_{42}','\theta_{52}')
     set(hleg1,'FontSize',30);

    
    subplot(1,3,2)
    set(gca,'fontsize',28)
    plot(T,Y(M*N*2+1:M*N*2+M,:),'-','LineWidth',1.5 )
    axis([0 end_time 0 400])
    xlabel('Time','FontSize',30)
    ylabel('s','FontSize',30)
    hleg2=legend('s_{1}','s_{2}')
    set(hleg2,'FontSize',30);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %here is a first version of welfare implementation
    %of course we can use another function to estimate welfare
    %welfare = 1 ./ (Y(M*N*2+1,:) + 1);
    %for i=2:M
    %    welfare = welfare + 1./ (Y(M*N*2+i,:) + 1);
    %end
    
    %with the above welfare, low stimuli drive the welfare to be high, 
    %while high stimuli don't have much importance
    %maybe the welfare should rather be determined by the high stimuli:
    sumstimuli = 0;
    for i=1:M
        sumstimuli = sumstimuli + Y(M*N*2+i,:);
    end 
    welfare = exp(-sumstimuli / (M*100));
    
    %subplot(3,2,4)
    %plot(T,welfare,'-')
    %axis([0 end_time 0 1])
    %title('welfare')
    
    % Linear colony growth in dependence of colony wealth
    N= N + welfare*beta;
    
    %subplot(3,2,5)
    %plot(T,N,'-')
    %axis([0 end_time 0 100])
    %title('N')
  
     subplot(1,3,3)
    set(gca,'fontsize',28)
    plot(T,welfare,'-','LineWidth',1.5)
    axis([0 end_time 0 1])
    xlabel('Time','FontSize',30)
    ylabel('W','FontSize',30)
    
 
 
    
    

end


function dy = ode_function(t,y)
    global M N phi zeta theta_min theta_max delta alpha p sigma
    dy = zeros(2*N*M+M,1);
    
    %theta matrix
    for j=1:M
        for i=1:N
            %factor = 0;
            %y((j-1)*N + i, 1);
            %if (y((j-1)*N + i, 1) > theta_min) && (y((j-1)*N + i, 1)<theta_max)
            %    factor = 1;
            %end 
            
            dy((j-1)*N + i, 1) = [(1 - y((j-1)*N + i +N*M, 1))*phi- y((j-1)*N + i +N*M, 1)*zeta];
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
            %gaussian_sp = sigma*randn;
            gaussian_sp = 0;
            
            dy(M*N+(j-1)*N + i, 1) = (y(M*N*2+j)^2/(y(M*N*2+j)^2+y((j-1)*N + i)^2)+0.0001)*(1-sum)-p*y(M*N+(j-1)*N + i, 1)+gaussian_sp;

           
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
    global M N theta_min theta_max sigma
    T = [];
    Y = [];
    T = [T,time(1)];
    Y = [Y,y0];
    
    
    y_eul = y0;
    t_eul = time(1);
    
    number = (time(2)-time(1))/h;
    
    for j=1:number
        %if(j==number/2)
            %theta matrix
         %   for k=1:M
          %      for i=1:N
           %         y_eul((k-1)*N + i) = 500;
            %        y_eul(N*M + (k-1)*N + i) = 0.1;
             %   end
            %end
        %end
        y_eul = y_eul + h*fun(t_eul,y_eul)+sqrt(h)*sigma*randn(size(y_eul));
        t_eul = t_eul + h ;
        j
        for k=1:M
            for i=1:N
                if(y_eul((k-1)*N + i)<theta_min)
                    y_eul((k-1)*N + i)=theta_min;
                end
                if(y_eul((k-1)*N + i)>theta_max)
                    y_eul((k-1)*N + i)=theta_max;
                end
                if(y_eul(N*M + (k-1)*N + i)< 0.0)
                    y_eul(N*M + (k-1)*N + i) = 0;
                end
                if(y_eul(N*M + (k-1)*N + i) > 1.0)
                    y_eul(N*M + (k-1)*N + i) = 1;
               end
            end
        
        for v=1:2*N*M+M
                if (y_eul(v)<=0)
                    y_eul(v)=0.0001;
                end 
        end

        end 
        
        Y = [Y, y_eul];
        T = [T, t_eul];
    end
    

end 