import numpy as np

sputter_time = 150  # Total time in ps during which the sputtering occurs
incr_d = 20  # Increase in radius in Å
incr_t = 10  # Time between radius increases in ps
incr_buffer = 20  # Time before the first radius increase in ps
radius = 70  # Initial radius in Å
area = np.pi*(radius/10)**2  # area of the sputtered area in nm^2
flux = 1  # Particle flux in particles/ps/nm^2
flux_area = flux*area    # flux per nm^2 in particles/ps


rng_seed = 214079644654894187569045748385860788528  # Random seed for the numpy generator, generated with SeedSequence.entropy
rng = np.random.default_rng()

current_time = 0
insert_times = []
fluxes = []
next_increase = incr_buffer

avg_fluxes = [[], [], [], [], []]  # Average fluxes for each radius increase


for j in range(1):
    while current_time < sputter_time:
        arrival_time = rng.exponential(scale=1/flux_area)
        current_time += arrival_time
        insert_times.append(current_time)

        # Dynamic flux - increase the radius of the sputtering region
        if current_time >= next_increase:
            fluxes.append([current_time, len(insert_times), area])
            radius = radius + incr_d
            area = np.pi*(radius/10)**2
            flux = flux - 0.075
            flux_area = flux*area
            next_increase += incr_t

    avg_flux = []
    for i in range(len(fluxes)):
        if i == 0:
            prev = [0, 0, 0]
        else:
            prev = fluxes[i-1]
        current = fluxes[i]

        t = current[0] - prev[0]
        p = current[1] - prev[1]
        a = current[2]

        if i <= 4:
            avg_fluxes[i].append(p/t/a)

for i in range(len(avg_fluxes)):
    print(f'Average flux: {np.mean(avg_fluxes[i])}')