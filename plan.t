
1) Estimation of parameters of models:
    a) N difference ()
    b) standard deviation ()
    c) seq depth (500-5000:500)
    d) number of objects (100, 500, 1000)
    c) number of replicates (5, 10, )

2) write scripts to generate sets
    set generation procedure:
        0) generate connected components based on poisson distribution
        1) get mean N for one species from each component based on normal distribution
        2) get StD for one species from each component based on normal distribution
        3) generate N
        4) generate rest of species from components
            a) multiply mean N by N_coefficient (0, 1]
            b) multiply StD by N_coefficient and V_coefficient(more then 1)
            c) generate N

    Total parameters for set generation:
        1) poisson lambda
        2) mean and StD for normal distributions of N and std of species
        3) seq depth
        4) number of species

    Parameters to disturb set:
       1) N of species with differential expression. (several values)
       2) Parameters of disturb (from uniform distribution)

3) generate sets based on parameter net
