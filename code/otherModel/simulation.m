%Main function
clearvars;

%Initialize the variables:
N=5; %Number of persons
M=3; %Number of tasks

taskvalue = ones(1,M);
chosenTask = zeros(N,1);
productivity = 5 + ones(N,M);% randn(N,M);
production = zeros(1,M);
money = zeros(N,1);
learning = 0.01;
forgetting = 0;

%Variables for simulation
timesteps=50;
stepSize=1;
output = zeros(N,1);

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

for t=0:timesteps
    %Update total production
    production = zeros(1,M);
    for i=1:N
        production(chosenTask(i,1)) =  production(chosenTask(i,1)) + productivity(i,chosenTask(i,1));
    end
    
    %Calculate who earns how much money
    for i=1:N
        money(i,1) = productivity(i,chosenTask(i,1)) / production(chosenTask(i,1));
    end
    
    %Update abilities
    for i=1:N
        for j=1:M
            if(chosenTask(i,1) == j)
                productivity(i,j) = min(10, productivity(i,j) + learning);
            else
                productivity(i,j) = max(1, productivity(i,j) - forgetting);
            end
        end
    end
    
    %Allow one person to change the job
    theChosenOne = mod(t,N) + 1;
    maxGain=0;
    bestTask=0;
    for j=1:M
        if(j == chosenTask(theChosenOne,1))
            gain = productivity(theChosenOne,j) / production(chosenTask(theChosenOne,1));
        else
            gain = productivity(theChosenOne,j) / (production(1,j)+productivity(theChosenOne,j));
        end
        if(gain > maxGain)
            maxGain=gain;
            bestTask=j;
        end
        [theChosenOne j gain]
    end
    chosenTask(theChosenOne,1) = bestTask;
    output = [output chosenTask];
end

