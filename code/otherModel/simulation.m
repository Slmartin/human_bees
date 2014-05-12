%Main function
clearvars;

%Set seed for random number generation
RandStream.setGlobalStream(RandStream('mt19937ar','seed',42))
RandStream.getGlobalStream()

%Initialize the variables:
N=7; %Number of persons
M=3; %Number of tasks

%Parameters of the model
abilitymu = 3;
abilitysd = 0.7;
prodmu = 2;
prodsd = 0.7;
boredommu=0.5;
boredomsd = 0.15;
learning = 0.01;
forgetting = 0.003;
boredomIncrease = 0.001;
boredomDecrease = 0.0003;
switchPossibilityFrequency = 0.003; %Meaning ~all 1/x days

%initialize random variables
abilities = abs(abilitymu + abilitysd*randn(N,1));
initialProductivity = abs(prodmu + prodsd*randn(N,M));
maximalBoredom = abs(boredommu + boredomsd*randn(N,M));


%Variables for simulation
simulationTime=10000;
dt = 1;

catastrophe = false;
catastropheTime = 3000;


%Other needed variables
taskvalue = ones(1,M);
chosenTask = ones(N,1);
productivity = initialProductivity;
boredom = zeros(N,M);
boredomAtChosenTask = zeros(N,1);
productivityAtChosenTask = zeros(N,1);
production = zeros(1,M);
money = zeros(N,1);
totalmoney = zeros(N,1);
taskTime = zeros(N,1);
moneyTime = zeros(N,1);
totalmoneyTime = zeros(N,1);
boredomTime = zeros(N,1);
productivityTime = zeros(N,1);
totalProduction = 0;
productionTime = 0;
alive = ones(N,1);

%Initialize dependent values
maximalProductivity = initialProductivity;
for i=1:N
    for j=1:M
        maximalProductivity(i,j) = maximalProductivity(i,j)*abilities(i,1);
    end
end
%maximalBoredom = initialProductivity*maxBoredomFactor;

%initialize choice: everyone choses task where productivity is best
maxProductivity = zeros(N,1);
for i=1:N
    for j=1:M
        if(productivity(i,j) > maxProductivity(i,1))
            chosenTask(i,1) = j;
            maxProductivity(i,1) = productivity(i,j);
        end
    end
end

timesteps = simulationTime/dt;
for t=0:timesteps
    %Update total production
    production = zeros(1,M);
    for i=1:N
        if(alive(i,1))
            production(chosenTask(i,1)) =  production(chosenTask(i,1)) + productivity(i,chosenTask(i,1));
        end
    end
    
    %Calculate who earns how much money and update total production
    totalProduction = 0;
    for i=1:N
        if(alive(i,1))
            money(i,1) = productivity(i,chosenTask(i,1)) / production(chosenTask(i,1));
            totalmoney(i,1) = totalmoney(i,1) + money(i,1);
            totalProduction = totalProduction + productivity(i,chosenTask(i,1));
        end
    end
    
    %Update abilities
    for i=1:N
        if(alive(i,1))
            for j=1:M
                if(chosenTask(i,1) == j)
                    productivity(i,j) = productivity(i,j) + (maximalProductivity(i,j)-productivity(i,j))*(learning*dt);
                    %productivity(i,j) = min(maximalProductivity(i,j), productivity(i,j) + learning*dt);
                    boredom(i,j) = boredom(i,j) + (maximalBoredom(i,j)-boredom(i,j))*(boredomIncrease*dt);
                    %boredom(i,j) = boredom(i,j) + (maximalBoredom(i,j)-boredom(i,j))*(boredomIncrease*boredomSensitivity(i,1)*dt);
                    %boredom(i,j) = min(maximalBoredom(i,j), boredom(i,j) + boredomIncrease*boredomSensitivity(i,1)*dt);
                else
                    productivity(i,j) = productivity(i,j) - (productivity(i,j)-initialProductivity(i,j))*(forgetting*dt);
                    %productivity(i,j) = max(initialProductivity(i,j), productivity(i,j) - forgetting*dt);
                    boredom(i,j) = boredom(i,j) - (boredom(i,j))*(boredomDecrease*dt);
                    %boredom(i,j) = boredom(i,j) - (boredom(i,j))*(boredomDecrease*boredomSensitivity(i,1)*dt);
                    %boredom(i,j) = max(0, boredom(i,j) - boredomDecrease*dt);
                end
            end
        end
    end
    
    %Allow one person to change the job
    for i=1:N
        if(alive(i,1))
            if(rand() < switchPossibilityFrequency*dt)
                maxGain=-1000;
                bestTask=0;
                for j=1:M
                    if(j == chosenTask(i,1))
                        gain = productivity(i,j) / production(chosenTask(i,1));
                    else
                        gain = productivity(i,j) / (production(1,j)+productivity(i,j));
                    end
                    gain = gain - boredom(i,j);
                    if(gain > maxGain)
                        maxGain=gain;
                        bestTask=j;
                    end
                end
                chosenTask(i,1) = bestTask;
            end
        end
    end
    
    %Update boredom
    for i=1:N
        if(alive(i,1))
            boredomAtChosenTask(i,1) = boredom(i,chosenTask(i,1));
            productivityAtChosenTask(i,1) = productivity(i,chosenTask(i,1));
        end
    end
    
    %Kill someone
    if(catastropheTime < t*dt)
        if catastrophe
            chosenTask(1,1) = -1;
            alive(1,1) = 0;
            totalmoney(1,1)=-1000;
            boredom(1,:)=-10;
            boredomAtChosenTask(1,1)=-10;
            productivityAtChosenTask(1,1)=-10;
            productivity(1,:)=-10;
            money(1,1)=-10;
        end
    end
    
    
    taskTime = [taskTime chosenTask];
    moneyTime = [moneyTime money];
    totalmoneyTime = [totalmoneyTime totalmoney];
    boredomTime = [boredomTime boredomAtChosenTask];
    productivityTime = [productivityTime productivityAtChosenTask];
    productionTime = [productionTime totalProduction];
end

time = 0:timesteps+1;

%For visibility, so that lines don't hide each other (superimposition); 
%TODO: find better way to make them all visible
for i=1:N
    for t=1:timesteps+1
        taskTime(i,t) = taskTime(i,t) - 0.1 + 0.2/N*i;
    end
end
fs=34;
th=5;
figure(1)
    set(gca,'fontsize',28)
    plot(time, taskTime, 'LineWidth',th)
    xlabel('Time','FontSize',fs)
    ylabel('Chosen Task number','FontSize',fs)
axis([0 timesteps 0.5 M+0.5])
set(gcf, 'PaperPosition', [0 0 28 14]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [28 14]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'taskno', 'pdf') %Save figure

figure(2)
    set(gca,'fontsize',28)
plot(time, moneyTime, 'LineWidth',th)
axis([0 timesteps 0 1.1])
    xlabel('Time','FontSize',fs)
    ylabel('Earned money','FontSize',fs)
set(gcf, 'PaperPosition', [0 0 28 14]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [28 14]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'money', 'pdf') %Save figure

figure(3)
    set(gca,'fontsize',28)
plot(time, productivityTime, 'LineWidth',th)
axis([0 timesteps 0 max(productivityTime(:))*1.1])
    xlabel('Time','FontSize',fs)
    ylabel('Productivity at chosen task','FontSize',fs)
    set(gcf, 'PaperPosition', [0 0 28 14]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [28 14]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'productivity', 'pdf') %Save figure

figure(4)
    set(gca,'fontsize',28)
plot(time, boredomTime, 'LineWidth',th)
axis([0 timesteps 0 max(boredomTime(:))+0.1])
    xlabel('Time','FontSize',fs)
    ylabel('Boredom at chosen task','FontSize',fs)
set(gcf, 'PaperPosition', [0 0 28 14]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [28 14]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'boredom', 'pdf') %Save figure

figure(5)
    set(gca,'fontsize',28)
plot(time, totalmoneyTime, 'LineWidth',th)
axis([0 timesteps 0 max(totalmoneyTime(:))*1.1])
    xlabel('Time','FontSize',fs)
    ylabel('Total amount of money earned so far','FontSize',fs)
set(gcf, 'PaperPosition', [0 0 28 14]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [28 14]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'totalmoney', 'pdf') %Save figure

figure(6)
subplot(4,2,2)
plot(time, productionTime)
title('Total production')


save()
%Idea for alternative representation (DIFFICULT!!)
%Have the tasks represented as locations on a 2D plot. Then each person
%(represented by a drawing of a person/worker or by a face)
%stands next to the task he is performing, and a pile of coins shows how
%much he earns. Then we see the workers moving to another point when they 
%change the task.