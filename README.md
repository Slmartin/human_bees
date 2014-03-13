# MATLAB Spring 2014 – Research Plan 

> * Group Name: human_bees
> * Group participants names: Schroeter Julien, Slusarczyk Martin, Vaucher Alain
> * Project Title: Human Bees

## General Introduction

Division of labor is an important aspect of society. It has been shown that even insect societies feature division of labor. It is therefore interesting to investigate the factors leading to this phenomenon. Moreover, a simple dynamic system can be used to gain insights which can be used to understand better the human society.

Theraulaz et al. have successfully developed a model describing the emergence of division of labor in insect societies. However their system is still very oversimplified. For instance they have studied systems exhibiting two tasks only. They do account that skills are developped in time, but they used fix population sizes, ignored aging, considered an isolated system without external disturbances and did not account for personal preferences.
It would therefore be interesting to explore further developments of their model.


## The Model

Our work is based on the model of Theraulaz et al., who investigated the division of labor in insect societies, based on variable response thresholds. 
In their model, the necessity to perform a task is described by a stimulus, to which each insect responds differently based on its respective threshold towards this task. These thresholds change in time depending on whether the insect is performing this task or not.
Additionally, on top of this model, we would like to introduce new factors describing aging and welfare. In an attempt to model a simplified human society. In our agent-based approach, aging will be modeled by a decrease in efficiency in time. Welfare is an observable that correlates negatively with the magnitude of the stimuli. Following the Penna Model, birth and death are introduced. The welfare influences the reproduction rate as well as the death probability. 
Further extensions could include the introduction of a currency, the "Honey Dollar", which would be connected to welfare and productivity. 

(Define dependent and independent variables you want to study. Say how you want to measure them.) (Why is your model a good abtraction of the problem you want to study?) (Are you capturing all the relevant aspects of the problem?)


## Fundamental Questions

We would like to study the dynamics of the model described above as well as the impact of the introduced factors. We would like to investigate the response of the model when exposed to external disturbances such as sudden death, wars, epidemics or vacation. 

- How does the system change when aging is introduced? Under which circumstances does the colony thrive respectively gets extinguished?
- How does the system behave upon temporary or definitive removal of a) random workers, b) young workers, c) workers fulfilling the same task?
- Is this model appropriate to explain a concept such as welfare?
- Can our model reproduce inegal wealth (money) distribution of the real world?


## Expected Results

Using our modified model, we still expect to reproduce the observation of labor division among bees. Upon application a disturbance, it is expected that the system reaches another equilibrium with a smaller welfare.
Without external constraints such as resource limit, it is expected that the population will increase indefinitely. 


## References 

Theraulaz, G., Bonabeau, E., Deneubourg, J.-L., Proc. R. Soc. Lond. B 265, 1998, 327-332.
Bonabeau, E., Theraulaz, G., Deneubourg, J.-L., Proc. R. Soc. Lond. B 263, 1996, 1565-1569.
Zoethout, K., Jager, W., Molleman, E., Simulation Modelling Practice and Theory 14, 2006, 342-359.
Penna, T.J.P., J. Stat. Phys. 78, 1995, 1629.


## Research Methods

Agent-Based Model


## Other

No dataset is going to be used.
