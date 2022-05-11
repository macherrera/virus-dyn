# Problem Set 3: Simulating the Spread of Disease and Virus Population Dynamics 

import random
import pylab
import numpy as np
from ps3b_precompiled_39 import *
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"

''' 
Begin helper code
'''

class NoChildException(Exception):
    """
    NoChildException is raised by the reproduce() method in the SimpleVirus
    and ResistantVirus classes to indicate that a virus particle does not
    reproduce. You can use NoChildException as is, you do not need to
    modify/add any code.
    """

'''
End helper code
'''

#
# PROBLEM 1
#
class SimpleVirus(object):

    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    def __init__(self, maxBirthProb, clearProb):
        """
        Initialize a SimpleVirus instance, saves all parameters as attributes
        of the instance.        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        clearProb: Maximum clearance probability (a float between 0-1).
        """
        self.maxBirthProb = maxBirthProb
        self.clearProb = clearProb


    def getMaxBirthProb(self):
        """
        Returns the max birth probability.
        """
        return getattr(self, "maxBirthProb")

    def getClearProb(self):
        """
        Returns the clear probability.
        """
        return getattr(self, "clearProb")

    def doesClear(self):
        """
        Stochastically determines whether this virus particle is cleared from the
        patient's body at a time step. 
        returns: True with probability self.getClearProb and otherwise returns
        False.
        """
        return random.random() < self.getClearProb()

    def reproduce(self, popDensity):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient and
        TreatedPatient classes. The virus particle reproduces with probability
        self.maxBirthProb * (1 - popDensity).
        
        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring SimpleVirus (which has the same
        maxBirthProb and clearProb values as its parent).         

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population.         
        
        returns: a new instance of the SimpleVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.               
        """
        new_virus_particle = SimpleVirus(self.getMaxBirthProb(),self.getClearProb())
        if random.random() < (self.getMaxBirthProb()*(1-popDensity)):
            return new_virus_particle
        else:
            raise NoChildException



class Patient(object):
    """
    Representation of a simplified patient. The patient does not take any drugs
    and his/her virus populations have no drug resistance.
    """    

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes.

        viruses: the list representing the virus population (a list of
        SimpleVirus instances)

        maxPop: the maximum virus population for this patient (an integer)
        """
        self.viruses = viruses
        self.maxPop = maxPop

    def getViruses(self):
        """
        Returns the viruses in this Patient.
        """
        return getattr(self, "viruses")

    def getMaxPop(self):
        """
        Returns the max population.
        """
        return getattr(self, "maxPop")

    def getTotalPop(self):
        """
        Gets the size of the current total virus population. 
        returns: The total virus population (an integer)
        """
        return len(self.viruses)

    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute the following steps in this order:
        
        - Determine whether each virus particle survives and updates the list
        of virus particles accordingly.   
        
        - The current population density is calculated. This population density
          value is used until the next call to update() 
        
        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.                    

        returns: The total virus population at the end of the update (an
        integer)
        """
        # First check which viruses survive, put in list surviving_virus
        surviving_virus = []
        for virus_particle in self.getViruses():
            if not virus_particle.doesClear():
                surviving_virus.append(virus_particle)
        new_pop_density = len(surviving_virus)/self.getMaxPop()

        # Then make a copy of surviving_virus for the next generation. Each s_virus_particle gets a chance to reproduce
        # If reproduction is successful, append the child. Set the viruses list to this next_gen_virus list, then
        # return the length of that list.
        next_gen_virus = surviving_virus[:]
        for s_virus_particle in surviving_virus:
            try:
                next_gen_virus.append(s_virus_particle.reproduce(new_pop_density))
            except NoChildException:
                pass
        setattr(self, "viruses", next_gen_virus)
        #self.viruses = next_gen_virus
        return len(next_gen_virus)

def test_virus(a,b):
    rona = SimpleVirus(a,b)
    mac = Patient([rona], 10000000000000)
    #print("viruses:", mac.getViruses())
    for n in range(100):
        print(mac.update())
        #print("viruses:", mac.getViruses())
    print("FINAL:", mac.getTotalPop())
    #print(mac.viruses)

#test_virus(0.4,0.4)


#
# PROBLEM 2
#
def simulationWithoutDrug(numViruses, maxPop, maxBirthProb, clearProb,
                          numTrials):
    """
    Run the simulation and plot the graph for problem 3 (no drugs are used,
    viruses do not have any drug resistance).    
    For each of numTrials trial, instantiates a patient, runs a simulation
    for 300 timesteps, and plots the average virus population size as a
    function of time.

    numViruses: number of SimpleVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: Maximum clearance probability (a float between 0-1)
    numTrials: number of simulation runs to execute (an integer)
    """
    # Make an initial list of numViruses, init_virus_pop.
    init_virus_pop = []
    for virus in range(numViruses):
        virus = SimpleVirus(maxBirthProb, clearProb)
        init_virus_pop.append(virus)
    # Make a list virus_populations to contain the lists for each patient. For each trial (patient) instantiate a
    # new_subject, with our init_virus_pop.
    virus_populations = []
    for trial in range(numTrials):
        new_subject = Patient(init_virus_pop, maxPop)
        new_virus_pop = []
        # Run 300 generations of the virus. Each time append the TotalPop to the new_virus_pop list. After 300
        # generations, append that list to the virus_populations list.
        for n in range(300):
            new_subject.update()
            new_virus_pop.append(new_subject.getTotalPop())
        #print("NEW PERSONS virus count:", new_virus_pop)
        virus_populations.append(new_virus_pop)
    #print(virus_populations)
    # To to get the average over generations: divide each value in the virus_populations list of lists by the
    # numTrials, put that in virus_pops_avgs. Then sum along the columns, to get the total_virus_avg, which is graphed.
    virus_pops_avgs = np.divide(virus_populations, numTrials)
    #print(virus_pops_avgs)
    total_virus_avg = [0] * 300

    for list in virus_pops_avgs:
        total_virus_avg += list
    x_val = [_ for _ in range(300)]
    #print(type(total_virus_avg.tolist()))
    pylab.plot(total_virus_avg.tolist(), label="SimpleVirus")
    pylab.title("SimpleVirus simulation")
    pylab.xlabel("Time Steps")
    pylab.ylabel("Average Virus Population")
    pylab.legend(loc="best")
    pylab.show()

#simulationWithoutDrug(1, 10, 1.0, 0.0, 1)


#
# PROBLEM 3
#
class ResistantVirus(SimpleVirus):
    """
    Representation of a virus which can have drug resistance.
    """   

    def __init__(self, maxBirthProb, clearProb, resistances, mutProb):
        """
        Initialize a ResistantVirus instance, saves all parameters as attributes
        of the instance.

        maxBirthProb: Maximum reproduction probability (a float between 0-1)       

        clearProb: Maximum clearance probability (a float between 0-1).

        resistances: A dictionary of drug names (strings) mapping to the state
        of this virus particle's resistance (either True or False) to each drug.
        e.g. {'guttagonol':False, 'srinol':False}, means that this virus
        particle is resistant to neither guttagonol nor srinol.

        mutProb: Mutation probability for this virus particle (a float). This is
        the probability of the offspring acquiring or losing resistance to a drug.
        """
        SimpleVirus.__init__(self, maxBirthProb, clearProb)
        self.resistances = resistances
        self.mutProb = mutProb

    def getResistances(self):
        """
        Returns the resistances for this virus.
        """
        return getattr(self, "resistances")

    def getMutProb(self):
        """
        Returns the mutation probability for this virus.
        """
        return getattr(self, "mutProb")

    def isResistantTo(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This method
        is called by getResistPop() in TreatedPatient to determine how many virus
        particles have resistance to a drug.       

        drug: The drug (a string)

        returns: True if this virus instance is resistant to the drug, False
        otherwise.
        """
        # Return the boolean resistance to a drug for this virus. If the drug is not in the dict, return False
        drug_dic = getattr(self, "resistances")
        #print(drug_dic, drug)
        try:
            return drug_dic[drug]
        except KeyError:
            return False


    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the TreatedPatient class.

        A virus particle will only reproduce if it is resistant to ALL the drugs
        in the activeDrugs list. For example, if there are 2 drugs in the
        activeDrugs list, and the virus particle is resistant to 1 or no drugs,
        then it will NOT reproduce.

        Hence, if the virus is resistant to all drugs
        in activeDrugs, then the virus reproduces with probability:      

        self.maxBirthProb * (1 - popDensity).                       

        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring ResistantVirus (which has the same
        maxBirthProb and clearProb values as its parent). The offspring virus
        will have the same maxBirthProb, clearProb, and mutProb as the parent.

        For each drug resistance trait of the virus (i.e. each key of
        self.resistances), the offspring has probability 1-mutProb of
        inheriting that resistance trait from the parent, and probability
        mutProb of switching that resistance trait in the offspring.       

        For example, if a virus particle is resistant to guttagonol but not
        srinol, and self.mutProb is 0.1, then there is a 10% chance that
        that the offspring will lose resistance to guttagonol and a 90%
        chance that the offspring will be resistant to guttagonol.
        There is also a 10% chance that the offspring will gain resistance to
        srinol and a 90% chance that the offspring will not be resistant to
        srinol.

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population       

        activeDrugs: a list of the drug names acting on this virus particle
        (a list of strings).

        returns: a new instance of the ResistantVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.
        """
        # Check for drug resistance. Default to fully_resistant, then if any drug is not resisted, change that
        # to False.
        tot_resistant = True
        for each_drug in activeDrugs:
            if not self.isResistantTo(each_drug):
                tot_resistant = False
                break
        # Virus reproduces if fully_resistant, AND stochastic test passes, else Exception.
        if tot_resistant and random.random() <= self.maxBirthProb * (1 - popDensity):
            # If reproduction is successful: make empty list of child_active_drugs, make virus_child with self properties
            # Look through each_drug in the Resistances list.
            child_active_drugs = []
            child_resistances = {}
            virus_child = ResistantVirus(self.maxBirthProb, self.clearProb, child_resistances, self.mutProb)
            for each_drug in self.getResistances():
                mutate = random.random() < self.getMutProb()
                parent_resistance = self.isResistantTo(each_drug)
                if mutate:
                    child_resistances[each_drug] = not parent_resistance
                else:
                    child_resistances[each_drug] = parent_resistance

                for key in child_resistances.keys():
                    if child_resistances[key]:
                        child_active_drugs.append(key)
            setattr(virus_child, "activeDrugs", child_active_drugs)
            setattr(virus_child, "resistances", child_resistances)
            return virus_child
        else:
            raise NoChildException


class TreatedPatient(Patient):
    """
    Representation of a patient. The patient is able to take drugs and his/her
    virus population can acquire resistance to the drugs he/she takes.
    """
    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered
        (which should initially include no drugs).              

        viruses: The list representing the virus population (a list of
        virus instances)

        maxPop: The  maximum virus population for this patient (an integer)
        """
        Patient.__init__(self, viruses, maxPop)
        self.drugs = []

    def addPrescription(self, newDrug):
        """
        Administer a drug to this patient. After a prescription is added, the
        drug acts on the virus population for all subsequent time steps. If the
        newDrug is already prescribed to this patient, the method has no effect.

        newDrug: The name of the drug to administer to the patient (a string).

        postcondition: The list of drugs being administered to a patient is updated
        """
        if newDrug not in self.drugs:
            self.drugs.append(newDrug)
        for particle in self.getViruses():
            particle.resistances[newDrug] = False


    def getPrescriptions(self):
        """
        Returns the drugs that are being administered to this patient.

        returns: The list of drug names (strings) being administered to this
        patient.
        """
        return getattr(self, "drugs")

    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed in
        drugResist.       

        drugResist: Which drug resistances to include in the population (a list
        of strings - e.g. ['guttagonol'] or ['guttagonol', 'srinol'])

        returns: The population of viruses (an integer) with resistances to all
        drugs in the drugResist list.
        """
        virus_particles = getattr(self, "viruses")

        num_resisting = 0
        for particle in virus_particles:
            virus_particle_resistance = True
            for drug in drugResist:
                #print("testing:", drug, "against:", particle.isResistantTo(drug))
                if not particle.isResistantTo(drug):
                    #print("setting virus_part_res to F:", particle.isResistantTo(drug))
                    virus_particle_resistance = False
            if virus_particle_resistance:
                #print("incrementing...")
                num_resisting += 1
        return num_resisting

    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:

        - Determine whether each virus particle survives and update the list of
          virus particles accordingly

        - The current population density is calculated. This population density
          value is used until the next call to update().

        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.
          The list of drugs being administered should be accounted for in the
          determination of whether each virus particle reproduces.

        returns: The total virus population at the end of the update (an
        integer)
        """
        # Make an empty list surviving_virus to contain the portion of current virus that remains. If there are current
        # virus particles, give each one a chance to clear and continue, append if it does. Calculate new_pop_density
        # with new population.
        surviving_virus = []
        for virus_particle in self.getViruses():
            #print(type(virus_particle))
            if not virus_particle.doesClear():
                surviving_virus.append(virus_particle)
        #print("len test:", len(surviving_virus), self.getMaxPop())
        new_pop_density = len(surviving_virus) / self.getMaxPop()

        next_gen_virus = surviving_virus[:]
        for s_virus_particle in surviving_virus:
            try:
                next_gen_virus.append(s_virus_particle.reproduce(new_pop_density, self.getPrescriptions()))
            except NoChildException:
                pass
        setattr(self, "viruses", next_gen_virus)
        return len(next_gen_virus)



def test_treated_patient():
    virus = ResistantVirus(1.0, 0.0, {None}, 0.0)
    patient = TreatedPatient([virus], 100)
    for i in range(99):
        print(patient.update())
    """virus1 = ResistantVirus(1.0, 0.0, {"drug1": True}, 0.0)
    virus2 = ResistantVirus(1.0, 0.0, {"drug1": False, "drug2": True}, 0.0)
    virus3 = ResistantVirus(1.0, 0.0, {"drug1": True, "drug2": True}, 0.0)
    patient = TreatedPatient([virus1, virus2, virus3], 100)
    print(patient.getResistPop(['drug1']))
    print(patient.getResistPop(['drug2']))
    print(patient.getResistPop(['drug1', 'drug2']))
    print("drug3:", patient.getResistPop(['drug3']))
    print(patient.getResistPop(['drug1', 'drug3']))
    print(patient.getResistPop(['drug1', 'drug2', 'drug3']))"""
#test_treated_patient()



#
# PROBLEM 4
#
def simulationWithDrug(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                       mutProb, numTrials):
    """
    Runs simulations and plots graphs for problem 5.

    For each of numTrials trials, instantiates a patient, runs a simulation for
    150 timesteps, adds guttagonol, and runs the simulation for an additional
    150 timesteps.  At the end plot the average virus population size
    (for both the total virus population and the guttagonol-resistant virus
    population) as a function of time.

    numViruses: number of ResistantVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: maximum clearance probability (a float between 0-1)
    resistances: a dictionary of drugs that each ResistantVirus is resistant to
                 (e.g., {'guttagonol': False})
    mutProb: mutation probability for each ResistantVirus particle
             (a float between 0-1). 
    numTrials: number of simulation runs to execute (an integer)
    
    """
    # Make an initial list of numViruses, init_virus_pop.
    init_virus_pop = []
    for virus in range(numViruses):
        virus = ResistantVirus(maxBirthProb, clearProb, resistances, mutProb)
        init_virus_pop.append(virus)
    # Make a list virus_populations to contain the lists for each patient. For each trial (patient) instantiate a
    # new_subject, with our init_virus_pop.
    virus_populations = []
    resistant_virus_populations = []
    for trial in range(numTrials):
        new_subject = TreatedPatient(init_virus_pop, maxPop)
        new_virus_pop = []
        new_resistant_virus_population = []
        # Run 300 generations of the virus. Each time append the TotalPop to the new_virus_pop list. After 300
        # generations, append that list to the virus_populations list.
        for n in range(300):
            if n == 150:
                new_subject.addPrescription("guttagonol")

            new_subject.update()

            new_virus_pop.append(new_subject.getTotalPop())
            #print(new_subject.getResistPop(["guttagonol"]))
            new_resistant_virus_population.append(new_subject.getResistPop(["guttagonol"]))
        # print("NEW PERSONS virus count:", new_virus_pop)
        virus_populations.append(new_virus_pop)
        resistant_virus_populations.append((new_resistant_virus_population))
    # print(virus_populations)

    # To to get the average over generations: divide each value in the virus_populations list of lists by the
    # numTrials, put that in virus_pops_avgs. Then sum along the columns, to get the total_virus_avg, which is graphed.
    virus_pops_avgs = np.divide(virus_populations, numTrials)
    resistance_virus_avgs = np.divide(resistant_virus_populations, numTrials)

    total_virus_avg = [0] * 300
    total_resistant_avg = [0] * 300

    for each_list in virus_pops_avgs:
        total_virus_avg += each_list
    for each_list in resistance_virus_avgs:
        #print(each_list)
        total_resistant_avg += each_list

    rounded_total_virus_avg = np.round(total_virus_avg, decimals=1)
    rounded_total_resistant_avg = np.round(total_resistant_avg, decimals=1)

    #print(rounded_total_virus_avg)
    #print(rounded_total_resistant_avg)
    pylab.plot(rounded_total_resistant_avg.tolist(), label="Resistant Virus")
    pylab.plot(rounded_total_virus_avg.tolist(), label="Total Virus")

    pylab.title("ResistantVirus simulation")
    pylab.xlabel("Time Steps")
    pylab.ylabel("Average Virus Population")
    pylab.legend(loc="best")
    pylab.show()


simulationWithDrug(1, 10, 1.0, 0.0, {}, 1.0, 5)

