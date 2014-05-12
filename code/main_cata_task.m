%---------------------------------
% Initialization of the variables
function recovery_time
global M N phi zeta theta_min theta_max delta alpha p sigma beta deads Somme Sommey task recovery_time count_bee
    
    VectorRecovery_task3=[];
    VectorRecovery_task2=[];
    for we=1:200
        %'start'
        N = 5;                      %Number of bees
        deads = 0;                  %Number of bees
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
        sigma=0.07;                 %sigma for gaussian distribution
        beta=50;                    %prefactor for growth rate

        dt=1;
        end_time = 10000;

        %dt=0.01;
        %end_time = 5000;


        theta_min=1;
        theta_max=1000;

        %[T,Y]=ode45(@ode_function, [0 3000], initial_condition);
        [T,Y]=euler_method(@ode_function, [0 end_time], initial_condition,dt);
        %var(Y(M*N*2+2,5001:6000))
        %var(Y(M*N*2+1,5001:6000))
        if (count_bee==2)
            VectorRecovery_task2 = [VectorRecovery_task2 round(recovery_time*dt)-end_time/2];
        else
            VectorRecovery_task3 = [VectorRecovery_task3 round(recovery_time*dt)-end_time/2];
        end
        we
    end
    
    
%fid = fopen('C:\Users\Julien\Desktop\Matlab_project\VectorRecovery_task3.txt','wt');  % Note the 'wt' for writing in text mode
%fprintf(fid,'%f\n',VectorRecovery_task3);  % The format string is applied to each element of a
%fclose(fid);
    
%fid = fopen('C:\Users\Julien\Desktop\Matlab_project\VectorRecovery_task2.txt','wt');  % Note the 'wt' for writing in text mode
%fprintf(fid,'%f\n',VectorRecovery_task2);  % The format string is applied to each element of a
%fclose(fid);    
    
end



function main_normal
    global M N phi zeta theta_min theta_max delta alpha p sigma beta deads Somme Sommey recovery_time count_bee

    Somme = [];
    Sommey =[];
    %'start';
    N = 5;                      %Number of bees
    deads = 0;                  %Number of bees
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
    sigma=0.07;                 %sigma for gaussian distribution
    beta=50;                    %prefactor for growth rate

    dt=1;
    end_time = 10000;

    %dt=0.01;
    %end_time = 5000;

    
    theta_min=1;
    theta_max=1000;
    
    %[T,Y]=ode45(@ode_function, [0 3000], initial_condition);
    [T,Y]=euler_method(@ode_function, [0 end_time], initial_condition,dt);
    

    %plot(T,Y(:,1:M*N),'-o')
    %plot(T,Y(:,13),'-r')
    subplot(6,2,1)
    plot(T,Y(1:M*N,:),'-')
    axis([0 end_time 0 1100])
    title('theta')
    
    %subplot(6,2,2)
    %plot(T,Y(M*N+1:M*N*2,:),'-')
    %axis([0 end_time -0.1 1.1])
    %title('x')
    
    subplot(6,2,2)
    plot(T,Y(M*N*2+1:M*N*2+M,:),'-')
    title('s')
    axis([0 end_time 0 10])
    
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
    
    %subplot(6,2,4)
    %plot(T,welfare,'-')
    %axis([0 end_time 0 1])
    %title('welfare')
    
    % Linear colony growth in dependence of colony wealth
    %N = N + welfare*beta;
    
    %subplot(3,2,5)
    %plot(T,N,'-')
    %axis([0 end_time 0 100])
    %title('N')
    
    %subplot(6,2,5)
    %plot(tsmovavg(Somme(1,:),'s',20,2),'-')
    %title('sum1')   
    

    
    %subplot(6,2,6)
    %plot(tsmovavg(Somme(2,:),'s',20,2),'-')
    %title('sum2')
    
 
    %subplot(6,2,7)
    %plot(Sommey(1,:),'-')
    %title('summande')   
    
    %subplot(6,2,8)
    %plot(Sommey(2,:),'-')
    %title('summande2')
    
    
    subplot(6,2,3)
    plot(T,Y(M*N+1,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x1')
    
    subplot(6,2,4)
    plot(T,Y(M*N+N+1,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x1_second')
    
    subplot(6,2,5)
    plot(T,Y(M*N+2,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x2')
    
    subplot(6,2,6)
    plot(T,Y(M*N+N+2,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x2_second')
    
    subplot(6,2,7)
    plot(T,Y(M*N+3,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x3')
    
    subplot(6,2,8)
    plot(T,Y(M*N+N+3,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x3_second')

    subplot(6,2,9)
    plot(T,Y(M*N+4,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x4')
    
    subplot(6,2,10)
    plot(T,Y(M*N+N+4,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x4_second')

    subplot(6,2,11)
    plot(T,Y(M*N+5,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x5')
    
    subplot(6,2,12)
    plot(T,Y(M*N+N+5,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x5_second')

end





function variance_data
%Seeds:
%mt19937ar random stream (current default)
%Seed: 0
%RandnAlg: Ziggurat

    %RandStream.getDefaultStream

    VectorVar2 = [];
    VectorVar32 = [];
    VectorMean2 = [];
    VectorMean32 = [];
  
    for j=1:500

        [vary, meany, i] = main_variance();
        j
        if(i==4)
            VectorVar32 = [VectorVar32 transpose(vary)];
            VectorMean32 = [VectorMean32 transpose(meany)];
        else
            VectorVar2 = [VectorVar2 transpose(vary)];
            VectorMean2 = [VectorMean2 transpose(meany)];
        end
        
    end
    

%fid = fopen('C:\Users\Julien\Desktop\Matlab_project\VectorVar32.txt','wt');  % Note the 'wt' for writing in text mode
%fprintf(fid,'%f\n',VectorVar32);  % The format string is applied to each element of a
%fclose(fid);

%fid = fopen('C:\Users\Julien\Desktop\Matlab_project\VectorVar2.txt','wt');  % Note the 'wt' for writing in text mode
%fprintf(fid,'%f\n',VectorVar2);  % The format string is applied to each element of a
%fclose(fid);

%fid = fopen('C:\Users\Julien\Desktop\Matlab_project\VectorMean32.txt','wt');  % Note the 'wt' for writing in text mode
%fprintf(fid,'%f\n',VectorMean32);  % The format string is applied to each element of a
%fclose(fid);

%fid = fopen('C:\Users\Julien\Desktop\Matlab_project\VectorMean2.txt','wt');  % Note the 'wt' for writing in text mode
%fprintf(fid,'%f\n',VectorMean2);  % The format string is applied to each element of a
%fclose(fid);
    

end









function [sol, meany, iy]=main_variance
global M N phi zeta theta_min theta_max delta alpha p sigma beta deads Somme Sommey task
    Somme = [];
    Sommey =[];
    %'start'
    N = 5;                      %Number of bees
    deads = 0;                  %Number of bees
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
    sigma=0.07;                 %sigma for gaussian distribution
    beta=50;                    %prefactor for growth rate

    dt=1;
    end_time = 10000;

    %dt=0.01;
    %end_time = 5000;

    
    theta_min=1;
    theta_max=1000;
    
    %[T,Y]=ode45(@ode_function, [0 3000], initial_condition);
    [T,Y]=euler_method(@ode_function, [0 end_time], initial_condition,dt);
 
    %subplot(3,1,1)
    %plot(T,Y(1:M*N,:),'-')
    %axis([0 end_time 0 1100])
    %title('theta')
    
    %subplot(3,1,2)
    %plot(T,Y(M*N+1:M*N*2,:),'-')
    %axis([0 end_time -0.1 1.1])
    %title('x')
    
    %subplot(3,1,3)    
    %plot(T,Y(M*N*2+1:M*N*2+M,:),'-')
    %title('s')
    %axis([0 end_time 0 10])
    

    i=0;
    elements2 = 1;
    elements3 = 2;
    %for j=1:N
    %a = tsmovavg(Y(M*N+j,3900:3999),'s',30,2);
    %    if(a(end)>0.2)
    %        i=i+1;
    %    end
    %end    
    %if (i==2)
    %    elements2 = 1;
    %    elements3 = 2;
    %end
    
    %iy=0;
    %for j=1:N*M
    %a = tsmovavg(Y(M*N+j,7900:7999),'s',30,2);
    %    if(a(end)>0.2)
    %        iy=iy+1;
    %    end
    %end
    
    
    if var(Y(M*N*2+elements3,5001:6000))>100 || var(Y(M*N*2+elements2,5001:6000))>100 
        iy=4;
        if(task==2) 
            elements2=2;
            elements3=1;
        end
        meany = [mean(Y(M*N*2+elements2,4500:4999)) mean(Y(M*N*2+elements3,4500:4999)) mean(Y(M*N*2+elements2,5001:6000)) mean(Y(M*N*2+elements3,5001:6000))];
        sol = [var(Y(M*N*2+elements2,4500:4999)) var(Y(M*N*2+elements3,4500:4999)) var(Y(M*N*2+elements2,5001:6000)) var(Y(M*N*2+elements3,5001:6000))];
    else
        iy=5;
        if (task==1)
            meany = [mean(Y(M*N*2+elements2,4500:4999)) mean(Y(M*N*2+elements3,4500:4999)) mean(Y(M*N*2+elements2,5001:6000)) mean(Y(M*N*2+elements3,5001:6000))];
            sol = [var(Y(M*N*2+elements2,4500:4999)) var(Y(M*N*2+elements3,4500:4999)) var(Y(M*N*2+elements2,5001:6000)) var(Y(M*N*2+elements3,5001:6000))];
        else
            elements2=2;
            elements3=1;
            
            meany = [mean(Y(M*N*2+elements2,4500:4999)) mean(Y(M*N*2+elements3,4500:4999)) mean(Y(M*N*2+elements2,5001:6000)) mean(Y(M*N*2+elements3,5001:6000))];
            sol = [var(Y(M*N*2+elements2,4500:4999)) var(Y(M*N*2+elements3,4500:4999)) var(Y(M*N*2+elements2,5001:6000)) var(Y(M*N*2+elements3,5001:6000))];
    
        end 
    end
    
    

   
end








function dy = ode_function(t,y)
    global M N phi zeta theta_min theta_max delta alpha p sigma deads Somme Sommey recovery_time
    dy = zeros(2*N*M+M,1);
    
    %theta matrix
    for j=1:M
        for i=1:N
            %factor = 0;
            %y((j-1)*N + i, 1);
            %if (y((j-1)*N + i, 1) > theta_min) && (y((j-1)*N + i, 1)<theta_max)
            %    factor = 1;
            %end 
            
            dy((j-1)*N + i, 1) = ((1 - y((j-1)*N + i +N*M, 1))*phi- y((j-1)*N + i +N*M, 1)*zeta);
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

           
            dy(M*N+(j-1)*N + i, 1) = (y(M*N*2+j)^2/(y(M*N*2+j)^2+y((j-1)*N + i)^2))*(1-sum)-p*y(M*N+(j-1)*N + i, 1)+gaussian_sp;
            
     

        end
    end 
    
    f=0;
    %s matrix
    for j=1:M
        %it makes no sense to have a negative value for the stiumuli
        sum=0;
     
        for i=1:N
            %it makes no sense to have a negative value for the stiumuli
            sum = sum + alpha/(N)*y(M*N+(j-1)*N+i);
        end
        
        if(f==0)
            f=1;
            ay=alpha/(N)*y(M*N+1);
            a=sum;
        else
            by=alpha/(N)*y(M*N+N+1);
            b=sum;
            
            cy=[ay by];
            c=[a b];
       
            
            Somme = [Somme, transpose(c)];
            Sommey=[Sommey,transpose(cy)];
        end
        
        
        if ~((y(2*M*N+j)<=0) && ((delta - sum)<0))
 
            dy(2*M*N+j)= dy(2*M*N+j) + delta - sum;
        end
    end
    
    
end


function [T,Y] = euler_method(fun,time, y0, h)
    global M N theta_min theta_max sigma deads task recovery_time count_bee
    T = [];
    Y = [];
    T = [T,time(1)];
    Y = [Y,y0];

    y_eul = y0;
    t_eul = time(1);
    
    number = (time(2)-time(1))/h;
    
    recovery_time = 0;
    for j=1:number
        
        %random breaks
        %if(round(1000*rand)==500)
        %    y_eul = cata_function3(y_eul,round((N-1)*rand)+1);
        %end
        
        %cata:take first one out
        %if(j==round((number/2)))
        %    task = 1;
        %    if (y_eul(M*N+1) < y_eul(M*N+1+N))
        %        task = 2;
        %    end
        %        
        %    y_eul = cata_function2(y_eul);
        %end
 
        %cata: take bees working actively on task 1 out
        %??difference in recovery: case count=2or3
        if(j==round((number/2)))
            count_bee =0;
            for bb=1:N
                if ((y_eul(N*M+bb)>0.2))
                    count_bee = count_bee+1;
                    y_eul = cata_function3(y_eul,bb);
                end
            end
                
            count_bee
        end
        
        
        
        if(j>round((number/2)+1000))
            if (y_eul(2*M*N+task)< 4)
               recovery_time = j;
               recovery=1;
               break
            end
        end
        
        
        y_eul = y_eul + h*fun(t_eul,y_eul)+sqrt(h)*sigma*randn(size(y_eul));
        t_eul = t_eul + h ;
       
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


function Y = cata_function(y0)
    global  deads N M
   
    %theta matrix
    deads = 2;
    i=1;
    Y=y0;
    for j=1:M
            Y((j-1)*N + i, 1) = 1000; 
    end
     
    %x matrix
    for j=1:M           
            Y(M*N+(j-1)*N + i, 1) = 0.000001;
    end
    
    
    i=2;
    Y=y0;
    for j=1:M
            Y((j-1)*N + i, 1) = 1000; 
    end
     
    %x matrix
    for j=1:M           
            Y(M*N+(j-1)*N + i, 1) = 0.000001;
    end
  
end


%restart all
function Y = cata_function2(y0)
    global  N M
    Y=y0;
    for k=1:M
        for i=1:N
            Y((k-1)*N + i) = 500;
            Y(N*M + (k-1)*N + i) = 0.1;
        end
    end

end 

%restart only one
function Y = cata_function3(y0,i)
    global  N M
    Y=y0;
    %i=1
    for k=1:M
            Y((k-1)*N + i) = 500;
            Y(N*M + (k-1)*N + i) = 0.1;
    end
    

    %i=2;
    for k=1:M
            Y((k-1)*N + i) = 500;
            Y(N*M + (k-1)*N + i) = 0.1;
    end

end
