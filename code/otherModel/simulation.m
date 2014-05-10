%Main function
clearvars;

%Initialize the variables:
N=7; %Number of persons
M=3; %Number of tasks

%Parameters of the model
abilities = ones(N,1) + abs(randn(N,1));
boredomSensitivity = abs(1+randn(N,1)/2.5);
initialProductivity = 1.0+abs(4*randn(N,M));
learning = 0.01;
forgetting = 0.003;
boredomIncrease = 0.001;
boredomDecrease = 0.0005;
maxBoredomFactor = 0.1;
switchPossibilityFrequency = 0.01; %Meaning ~all 1/x days

%Variables for simulation
simulationTime=10000;
dt = 1;

%Other needed variables
taskvalue = ones(1,M);
chosenTask = ones(N,1);
productivity = initialProductivity;
boredom = zeros(N,M);
boredomAtChosenTask = zeros(N,1);
production = zeros(1,M);
money = zeros(N,1);
taskTime = zeros(N,1);
moneyTime = zeros(N,1);
boredomTime = zeros(N,1);
totalProduction = 0;
productionTime = 0;

%Initialize dependent values
maximalProductivity = initialProductivity;
for i=1:N
    for j=1:M
        maximalProductivity(i,j) = maximalProductivity(i,j)*abilities(i,1);
    end
end
maximalBoredom = initialProductivity*maxBoredomFactor;

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
        production(chosenTask(i,1)) =  production(chosenTask(i,1)) + productivity(i,chosenTask(i,1));
    end
    
    %Calculate who earns how much money and update total production
    totalProduction = 0;
    for i=1:N
        money(i,1) = productivity(i,chosenTask(i,1)) / production(chosenTask(i,1));
        totalProduction = totalProduction + productivity(i,chosenTask(i,1));
    end
    
    %Update abilities
    for i=1:N
        for j=1:M
            if(chosenTask(i,1) == j)
                productivity(i,j) = min(maximalProductivity(i,j), productivity(i,j) + learning*dt);
                boredom(i,j) = boredom(i,j) + (maximalBoredom(i,j)-boredom(i,j))*(boredomIncrease*boredomSensitivity(i,1)*dt);
                %boredom(i,j) = min(maximalBoredom(i,j), boredom(i,j) + boredomIncrease*boredomSensitivity(i,1)*dt);
            else
                productivity(i,j) = max(initialProductivity(i,j), productivity(i,j) - forgetting*dt);
                boredom(i,j) = boredom(i,j) - (boredom(i,j))*(boredomDecrease*boredomSensitivity(i,1)*dt);
                %boredom(i,j) = max(0, boredom(i,j) - boredomDecrease*dt);
            end
        end
    end
    
    %Allow one person to change the job
    for i=1:N
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
    
    %Update boredom
    for i=1:N
        boredomAtChosenTask(i,1) = boredom(i,chosenTask(i,1));
    end
    
    taskTime = [taskTime chosenTask];
    moneyTime = [moneyTime money];
    boredomTime = [boredomTime boredomAtChosenTask];
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

subplot(3,2,1)
plot(time, taskTime)
axis([0 timesteps 0 M+1])
title('Chosen task number');

subplot(3,2,3)
plot(time, moneyTime)
axis([0 timesteps 0 1.1])
title('Money earned');

subplot(3,2,5)
plot(time, boredomTime)
axis([0 timesteps 0 max(boredomTime(:))+0.1])
title('Boredom at chosen task');

subplot(3,2,2)
plot(time, productionTime)
title('Total production')
%Idea for alternative representation (DIFFICULT!!)
%Have the tasks represented as locations on a 2D plot. Then each person
%(represented by a drawing of a person/worker or by a face)
%stands next to the task he is performing, and a pile of coins shows how
%much he earns. Then we see the workers moving to another point when they 
%change the task.
