###Measurement data of a Predator-Prey model
### Column 1: population size of prey
### Column 2: population size of predators
function predatorpreydataset1()
    #create the data
    timepoints = linspace(0,50,6)
    measurements =
        [51.7002527839933	3.62886163423551;
         88.4059566464213	11.8592208786134;
         50.5971643551162	33.1299499686354;
         25.2455746232269	30.1443180597855;
         42.2660932030265	9.00691042344652;
         51.4088102060744	29.8816967299902]

    #return the data
    timepoints,measurements
end
