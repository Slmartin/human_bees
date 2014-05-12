
%---------------------------------
% Initialization of the variables


function main
    global M N phi zeta theta_min theta_max delta alpha p sigma beta Y_final N_saved N_total
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
    zeta=4.5;                     %learning %4
    phi=1;                      %forgetting
    sigma=0.07;                  %sigma for gaussian distribution
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
    
  


    
    subplot(9,2,1)
    %plot(T,Y(1:M*N,:),'-')
    plot(T,Y_final(1:M*N_saved,:))
    %axis([0 end_time 0 1100])
    title('theta')
 
    subplot(9,2,2)
    plot(T(1:end_time),N_total,'-')
    axis([0 end_time 0 70])
    title('N')
    
    subplot(9,2,3)
    plot(T,Y(M*N*2+1:M*N*2+M,:),'-')
    title('s')
    axis([0 end_time 0 400])
    
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
    
    subplot(9,2,4)
    plot(T,welfare,'-')
    axis([0 end_time 0 1])
    title('welfare')
    
    
    subplot(9,2,5)
    plot(T,Y_final(M*N_saved+1,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x1')
    
    subplot(9,2,6)
    plot(T,Y_final(M*N_saved+2,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x1_second')
    
    subplot(9,2,7)
    plot(T,Y_final(M*N_saved+3,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x2')
    
    subplot(9,2,8)
    plot(T,Y_final(M*N_saved+4,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x2_second')
    
    subplot(9,2,9)
    plot(T,Y_final(M*N_saved+5,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x3')
    
    subplot(9,2,10)
    plot(T,Y_final(M*N_saved+6,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x3_second')
    
    
    subplot(9,2,11)
    plot(T,Y_final(M*N_saved+7,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x4')
    
    subplot(9,2,12)
    plot(T,Y_final(M*N_saved+8,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x4_second')
    
    
    subplot(9,2,13)
    plot(T,Y_final(M*N_saved+9,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x5')
    
    subplot(9,2,14)
    plot(T,Y_final(M*N_saved+10,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x5_second')
    
    
    subplot(9,2,15)
    plot(T,Y_final(M*N_saved+11,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x6')
    
    subplot(9,2,16)
    plot(T,Y_final(M*N_saved+12,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x6_second')
    
    
    subplot(9,2,17)
    plot(T,Y_final(M*N_saved+11,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x7')
    
    subplot(9,2,18)
    plot(T,Y_final(M*N_saved+12,:),'-')
    axis([0 end_time -0.1 1.1])
    title('x7_second')
    

end


function dy = ode_function(t,y)
    global M N phi zeta theta_min theta_max delta alpha p sigma N_saved N_total
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
    global M N theta_min theta_max sigma N_saved Y_final N_total
    T = [];
    Y = [];
    N_total = [];
    T = [T,time(1)];
    Y = [Y,y0];
    
    Y_final=[];
    
    
    y_eul = y0;
    t_eul = time(1);
    
    number = (time(2)-time(1))/h;
    N_saved =0;
    for j=1:number
        
    factor=1.1;
            
    sumstimuli = 0;
    for i=1:M
        sumstimuli = sumstimuli + Y(M*N*2+i,:);
    end 
    welfare = exp(-sumstimuli / (M*100));
    
        ff = 1000*rand;
        %%kill someone
        if(ff>950) 
            tt = rand*factor;
            if(welfare>tt)
             [Y,y_eul]=add_bee(Y,y_eul);
            end

        end
        
        ff2 = 1000*rand;
        if(ff>950) 
           tt = rand*factor;
           if(welfare<tt)
                kill_number = round(rand*N);
                [Y,y_eul]=kill_bee(Y,y_eul,1);
           end
        end
        

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
        N_total = [N_total, N];

    end
    
    aaa=size(Y_final);
    aa=size(Y);
    
          if (aaa(1)>0)
              final_plus = [Y_final, zeros(aaa(1), aa(2)-aaa(2))];
          else
              final_plus = Y_final;
          end

           fin = final_plus(1:N_saved*M,:) ;
           
           for jj=1:N
           for k=1:M
                fin = [fin;Y(jj+(k-1)*N,:)];
           end
           end
           
           %%x
           fin = [fin;final_plus(M*N_saved+1:2*N_saved*M,:) ];
           for jj=1:N
           for k=1:M
                fin = [fin;Y(N*M+jj+(k-1)*N,:)];
           end
           end
           
            Y_final = fin;
            N_saved=N_saved+N;
 
  
end 



function [Y_sol,y_eul_sol] = kill_bee(Y,y_eul,kill_number)
    global M N theta_min theta_max sigma N_saved Y_final
            
        aaa=size(Y_final);
        aa=size(Y);

          %%theta
          if (aaa(1)>0)
              final_plus = [Y_final, zeros(aaa(1), aa(2)-aaa(2))];
          else
              final_plus = Y_final;
          end

           fin = final_plus(1:N_saved*M,:) ;
          
           for k=1:M
                fin = [fin;Y(kill_number+(k-1)*N,:)];
           end
           
           %%x
           fin = [fin;final_plus(M*N_saved+1:2*N_saved*M,:) ];
           for k=1:M
                fin = [fin;Y(N*M+kill_number+(k-1)*N,:)];
           end
           
            Y_final = fin;
 
           for k=1:M
                [Y,ps] =  removerows(Y,'ind',N*M+kill_number+(M-k)*N);
                [y_eul,ps] =  removerows(y_eul,'ind',N*M+kill_number+(M-k)*N);
           end
           for k=1:M
                [Y,ps] =  removerows(Y,'ind',kill_number+(M-k)*N);
                [y_eul,ps] =  removerows(y_eul,'ind',kill_number+(M-k)*N);
           end
           
            N_saved  = N_saved +1;
            N=N-1;  
            
            Y_sol = Y;
            y_eul_sol = y_eul;
    
end


function [Y_sol,y_eul_sol] = add_bee(Y,y_eul)
    global M N theta_min theta_max sigma N_saved Y_final
        
 
    aa=size(Y);

          %%theta

            Y_sol = [];
            y_eul_sol = [];
          
           for k=1:M
                Y_sol = [Y_sol; Y(1+N*(k-1):N*(k),:)] ;
                y_eul_sol = [y_eul_sol; y_eul(1+N*(k-1):N*(k),:)];

                Y_sol = [Y_sol; zeros(1,aa(2))] ;
                y_eul_sol = [y_eul_sol; 500];
           end
           
           %%x
           a = size(Y_sol)
           b = size(y_eul_sol)

           for k=1:M
                Y_sol = [Y_sol; Y(N*M+N*(k-1)+1:N*M+N*(k),:)] ;
                y_eul_sol = [y_eul_sol; y_eul(N*M+N*(k-1)+1:N*M+N*(k),:)];

                Y_sol = [Y_sol; zeros(1,aa(2))] ;
                y_eul_sol = [y_eul_sol; 0];
                
           end
           
                Y_sol = [Y_sol; Y(2*N*M+1:2*N*M+2,:)] ;
                y_eul_sol = [y_eul_sol; y_eul(2*N*M+1:2*N*M+2,:)];      
    
           N = N + 1;
end
