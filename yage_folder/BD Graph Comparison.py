#!/usr/bin/env python
# coding: utf-8

# In[37]:


import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def ode_system(t, M, k1, k2):
    return -k1 * M + k2

# Define the initial conditions and parameters
initial_condition = 0.0
k1 = 0.1
k2 = 0.8

# Set the integration time span
time_span = (0.0, 100.0)

# Call solve_ivp to solve the ODE
solution = solve_ivp(ode_system, time_span, [initial_condition], args=(k1, k2))

# Extract the solution
time_points = solution.t
solution_values = solution.y[0]

# Plot the solution
plt.plot(time_points, solution_values,'--')
plt.xlabel('Time')
plt.ylabel('M')
plt.title('Solution of dM/dt = -k1 * M + k2')
plt.show()


# In[49]:


import numpy as np
import matplotlib.pyplot as plt

# Define birth and death rates
birth_rate = 0.8
death_rate = 0.1

# Define initial population size
initial_population = 0

# Define simulation parameters
t_max = 100
num_trajectories = 10

# Function to perform Gillespie simulation
def gillespie_simulation(birth_rate, death_rate, initial_population, t_max):
    time_points = [0]
    population_sizes = [initial_population]
    
    while time_points[-1] < t_max:
        # Calculate propensities
        prop_birth = birth_rate
        prop_death = death_rate * population_sizes[-1]
        total_propensity = prop_birth + prop_death
        
        # Generate random numbers
        r1, r2 = np.random.random(size=2)
        
        # Calculate time increment
        tau = -np.log(r1) / total_propensity
        
        # Determine the reaction
        if r2 < prop_birth / total_propensity:
            population_sizes.append(population_sizes[-1] + 1)  # Birth event
        else:
            population_sizes.append(population_sizes[-1] - 1)  # Death event
        
        # Update time
        time_points.append(time_points[-1] + tau)
    
    return time_points, population_sizes

# Perform Gillespie simulation for multiple trajectories
trajectories = []
for _ in range(num_trajectories):
    time_points, population_sizes = gillespie_simulation(birth_rate, death_rate, initial_population, t_max)
    trajectories.append((time_points, population_sizes))
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# Define the function representing the ODE
def model(M, t, k1, k2):
    dMdt = -k1 * M + k2
    return dMdt

# Define the parameters
k1 = 0.1
k2 = 0.8

# Define the time points
t = np.linspace(0, 100, 100)

# Solve the ODE
M = odeint(model, 0, t, args=(k1, k2))

# Plot the trajectories



plt.figure(figsize=(8, 6))
for i, (time_points, population_sizes) in enumerate(trajectories):
    plt.step(time_points, population_sizes, where='post', label=f'Trajectory {i+1}', linewidth=0.5)
plt.xlabel('Time[sec]')
plt.ylabel('Population size')
#plt.title('Birth and Death Process - Gillespie Simulation')

mean, = plt.plot(t, M,'--',c= 'black',label = 'mean')
plt.legend(handles=[mean], fontsize = '20')


# In[ ]:





# In[ ]:





# In[ ]:





# In[36]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
k = np.arange(0, 25)

#print Poission distribution
pmf = poisson.pmf(k, mu=8)
pmf = np.round(pmf, 5)

birth_rate = 0.8
death_rate = 0.1

# Define initial population size
initial_population = 0
# Define simulation parameters
t_max = 100
def gillespie_simulation(birth_rate, death_rate, initial_population, t_max):
    time_points = [0]
    population_sizes = [initial_population]
    
    while time_points[-1] < t_max:
        # Calculate propensities
        prop_birth = birth_rate
        prop_death = death_rate * population_sizes[-1]
        total_propensity = prop_birth + prop_death
        
        # Generate random numbers
        r1, r2 = np.random.random(size=2)
        
        # Calculate time increment
        tau = -np.log(r1) / total_propensity
        
        # Determine the reaction
        if r2 < prop_birth / total_propensity:
            population_sizes.append(population_sizes[-1] + 1)  # Birth event
        else:
            population_sizes.append(population_sizes[-1] - 1)  # Death event
        
        # Update time
        time_points.append(time_points[-1] + tau)
   # print(population_sizes)
    return time_points, population_sizes
    
def call_runtime_against_k():
    iters = 100000
    #fig, ax = plt.subplots(figsize =(10, 7))
    end_vals = []
    for i in range(iters):
        e, X = gillespie_simulation(birth_rate, death_rate, initial_population, t_max)
        end_vals.append(X[-1])
    max_val = int(np.max(end_vals))
    bins = np.array([i for i in range(max_val + 2)]) - 0.5
    kwargs = dict(density=True, linewidth=1.5, bins=bins, ec='black')
    #print(e)
    # p_points = np.linspace(-0.5, max_val + 0.5, 100)
    plt.xlabel("Molecule number")
    plt.ylabel("Stationary distribition")
    #ax.set_title("Stationary distribution")
    plt.hist(end_vals,
                  #density=True,
                  range=(0, max_val + 1),
                  fc='white',
                  **kwargs,label='Gillespie \n SSA')
    #mu_1 = DB_X0[0] * np.exp(- DB_rates[0] * 5)
    #mu_2 = sum([a * b for a, b in enumerate(pdf[0])])
   #print(end_vals) 
    #plt.show()
call_runtime_against_k()
plt.plot(k, pmf, 'r-', label='master\n equation')
plt.legend(fontsize=20)
#plt.rcParams.update({'font.size':25})
plt.show()


# In[ ]:





# In[ ]:




