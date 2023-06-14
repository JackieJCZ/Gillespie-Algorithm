#!/usr/bin/env python
# coding: utf-8

# In[53]:


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


#solving recurrence equation 
import numpy as np
import matplotlib.pyplot as plt

# Recurrence formula
def recurrence_formula(n, phi_prev, phi_curr, phi_prev2, k1, k2):
    return (k1 * n * phi_curr + k2 * phi_curr - k2 * phi_prev2) / (k1 * (n + 1))

# Parameters
k1 = 0.1
k2 = 0.8
n_max = 22  # Maximum value of n for the plot

# Initial conditions
phi_0 = 1
phi_1 = (k2 / k1) * phi_0

# Arrays to store values
n_values = np.arange(n_max + 1)
phi_values = np.zeros(n_max + 1)

# Compute phi values using the recurrence formula
phi_values[0] = phi_0
phi_values[1] = phi_1

for n in range(1, n_max):
    phi_values[n + 1] = recurrence_formula(n, phi_values[n], phi_values[n], phi_values[n - 1], k1, k2)

# Compute the sum of phi values
phi_sum = np.sum(phi_values)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(n_values, phi_values/phi_sum , 'r-', label='Numerical Solution')

# Add vertical line at the sum of phi(n)
##plt.axvline(x=phi_sum, color='k', linestyle='--', label='Sum of phi(n)')

plt.xlabel('number of molecules')
plt.ylabel('stationary distribution')
plt.title('Numerical Solution of the Recurrence Formula')
plt.grid(True)
plt.legend()
plt.show()


# In[48]:


import numpy as np
import matplotlib.pyplot as plt

# Recurrence formula
def recurrence_formula(n, phi_prev, phi_curr, phi_prev2, k1, k2):
    return (k1 * n * phi_curr + k2 * phi_curr - k2 * phi_prev2) / (k1 * (n + 1))

# Parameters
k1 = 0.1
k2 = 0.8
n_max = 25  # Maximum value of n for the plot

# Initial conditions
phi_0 = 1
phi_1 = (k2 / k1) * phi_0

# Arrays to store values
n_values = np.arange(n_max + 1)
phi_values = np.zeros(n_max + 1)

# Compute phi values using the recurrence formula
phi_values[0] = phi_0
phi_values[1] = phi_1

for n in range(1, n_max):
    phi_values[n + 1] = recurrence_formula(n, phi_values[n], phi_values[n], phi_values[n - 1], k1, k2)

# Compute the sum of phi values
phi_sum = np.sum(phi_values)

# Define birth and death rates
birth_rate = 0.8
death_rate = 0.1

# Define initial population size
initial_population = 0

# Define simulation parameters
t_max = 10000
num_repeats = 100

# Define recording interval
record_interval = 1  # Record population size every second

# Function to perform Gillespie simulation
def gillespie_simulation(birth_rate, death_rate, initial_population, t_max, record_interval):
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
    
    # Record population size at specified intervals
    time_points_recorded = []
    population_sizes_recorded = []
    for i in range(0, len(time_points), record_interval):
        time_points_recorded.append(time_points[i])
        population_sizes_recorded.append(population_sizes[i])
    
    return time_points_recorded, population_sizes_recorded

# Perform Gillespie simulation multiple times
max_population = 0
recorded_population_sizes = []
for _ in range(num_repeats):
    time_points, population_sizes = gillespie_simulation(birth_rate, death_rate, initial_population, t_max, record_interval)
    max_population = max(max_population, max(population_sizes))
    recorded_population_sizes.extend(population_sizes)

# Compute the maximum population size
max_population = max(recorded_population_sizes)

# Compute the histogram
population_counts, _ = np.histogram(recorded_population_sizes, bins=np.arange(max_population + 2))

# Normalize the histogram
phi_n_simulated = population_counts / len(recorded_population_sizes)

plt.figure(figsize=(10, 6))
plt.plot(n_values, phi_values / phi_sum, 'r-', label='master equation')
plt.bar(np.arange(max_population + 1), phi_n_simulated, alpha=0.5, label='Gillespie Simulation')

#adding the dashed vertical line at x=8
x_coord = np.where(n_values == 8)[0][0]
plt.axvline(x=x_coord, color='k', linestyle='--', dashes=(5, 5))

plt.xlabel('Number of molecules (n)')
plt.ylabel('Stationary distribution (phi(n))')
#plt.title('Comparison of Numerical Solution and Gillespie Simulation')

plt.legend(fontsize = '20')
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




