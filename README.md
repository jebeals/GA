# GA
Genetic Algorithims for Global Pseduo-Optimization

This is a Genetic algorithim script in its infancy that tests the algorithim over Ackley's function. Chromosomes were combined to be one large string
of binary, but it still work with Random crossover and mutation; improvements can and will be made to seperate chromosomes. One other thing worthy to metnion
is that the code only manipuately integer binary WRT MATLAB's dec2bin function. This also means that only positive values make sense too. For this reason, 
an alterior sign logic loop had to be included, and if one were to impleement this algorithim over an asymmetric domain about the origin, the architecture would
have to be re-tought. That being said, there is a way around the dec2bin only affecting integers (via expanding and contracting the decimals by factors of 10)
but the negative domain still reamins an issue.

A pretty good script to familiarze oneself with the high-level idea of GA in general.
