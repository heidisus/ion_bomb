import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

filename = 'tables/explin_300K_f25_900ps_100eV_r70A_incr5A.tsv'
# filename = 'tables/explin_300K_f26_150ps_100eV_r70A_incr20A.tsv'
# title = 'Temperature'
# ylabel = 'T (K)'
# xlabel = 'Time (ps)'

# Determine the # of lines after output data, then read the output
rows_skipped = 0

info = []
# Find number of lines to skip at the end of the file
with open(filename, 'r') as input:
    lines = input.readlines()
    for line_idx in range(len(lines)-1,0, -1):
        try:
            int(lines[line_idx].split()[0])
            break
        except:
            rows_skipped += 1
    info = lines[len(lines)-rows_skipped:]  # Add the data printed at the end of the run to the info list

print(rows_skipped)
df = pd.read_csv(filename, delim_whitespace=True, engine='python', skipfooter=rows_skipped) 

# df['norm_atoms'] = df['Atoms']-df['Atoms'][0]


# fig, ax = plt.subplots()
# df.plot(x='Time', y=['Temp', 'c_tbox'], label=['Entire box', 'Initial surface region'])
# df.plot(x='Time', y=['Atoms'])


# plt.title(title)
# plt.xlabel(xlabel)
# plt.ylabel(ylabel)
# plt.grid()

# plt.legend()
# plt.show()

# Time, flux, area, radius

# f25:
flux_theoretical = np.array([[0.0, 70.0], 
                             [150.0,70.0], 
                             [150.09256233564133,75.0], 
                             [200.0230692044839, 80.0], 
                             [250.06529160135344,85.0], 
                             [300.0212840814889, 90.0], 
                             [350.0098710813766, 95.0], 
                             [400.070339063182, 100.0], 
                             [450.06288423099596,105.0 ], 
                             [500.08622445330286,110.0 ], 
                             [550.0628783089085, 115.0], 
                             [600.029072242221, 120.0], 
                             [650.0432175523679, 125.0], 
                             [700.0432141359423, 130.0], 
                             [750.0355717971834, 135.0], 
                             [800.2756883138152, 140.0], 
                             [850.0253312700336, 145.0], 
                             [900.2600055313004, 150.0 ], ])

flux_effective =  [[150.09256233564133, 2300, 153.93804002589985], 
                   [200.0230692044839, 3129, 176.71458676442586], 
                   [250.06529160135344, 4017, 201.06192982974676], 
                   [300.0212840814889, 4923, 226.98006922186255], 
                   [350.0098710813766, 5901, 254.46900494077323], 
                   [400.070339063182, 6858, 283.5287369864788], 
                   [450.06288423099596, 7805, 314.1592653589793], 
                   [500.08622445330286, 8769, 346.3605900582747], 
                   [550.0628783089085, 9721, 380.132711084365], 
                   [600.029072242221, 10595, 415.4756284372501], 
                   [650.0432175523679, 11397, 452.3893421169302], 
                   [700.0432141359423, 12074, 490.8738521234052], 
                   [750.0355717971834, 12694, 530.929158456675], 
                   [800.2756883138152, 13180, 572.5552611167398], 
                   [850.0253312700336, 13486, 615.7521601035994], 
                   [900.2600055313004, 13605, 660.519855417254],]

# f26
# flux_theoretical = np.array([[0.0, 70.0],
#                              [20.0, 70.0],
#                              [20.008035791656482, 90.0],
#                              [30.006000726769038, 110.0],
#                              [40.00032266854119, 130.0],
#                              [50.00281293613858, 150.0],
#                              [60.000435042199726, 170.0],
#                              [70.00050208713506, 190.0],
#                              [80.00139406374468, 210.0],
#                              [90.00130363244152, 230.0],
#                              [100.00003210511379, 250.0],
#                              [110.00046491683524, 270.0],
#                              [120.00288831995648, 290.0],
#                              [130.00103368603243, 310.0],
#                              [140.00036459115606, 330.0],
#                              [150.00018770772644, 350.0]])


# flux_effective = [[20.008035791656482, 3095, 153.93804002589985], [30.006000726769038, 5462, 254.46900494077323], [40.00032266854119, 8737, 380.132711084365], [50.00281293613858, 12910, 530.929158456675], [60.000435042199726, 18074, 706.8583470577034], [70.00050208713506, 24094, 907.9202768874502], [80.00139406374468, 30787, 1134.1149479459152], [90.00130363244152, 38016, 1385.4423602330987], [100.00003210511379, 45524, 1661.9025137490005], [110.00046491683524, 53096, 1963.4954084936207], [120.00288831995648, 60429, 2290.221044466959], [130.00103368603243, 66988, 2642.079421669016], [140.00036459115606, 72462, 3019.0705400997913], [150.00018770772644, 76426, 3421.194399759285]]

# 1.0048708428926627 0.09954574735587428
def flux_effective_func(fluxes):
    t_result = [0.0]
    f_result = [0.09954574735587428*10**26]
    for i in range(len(fluxes)):
        if i == 0:
            prev = [0, 0, 0]
        else:
            prev = fluxes[i-1]
        current = fluxes[i]

        t = current[0] - prev[0]
        p = current[1] - prev[1]
        a = current[2]
        try:
            t_result.append(current[0])
            f_result.append((p/t/a)*10**26)
        except ZeroDivisionError:
            continue
    # print(f'Average flux: {np.mean(f_result)}')

    # for i in range(len(f_result)):
    #     print(f'{t_result[i]} {f_result[i]}')

    return t_result, f_result, np.mean(f_result)


rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
plt.rcParams.update({'font.size': 11})

fig, [ax, ax2] = plt.subplots(2, 1)
ax.plot(flux_effective_func(flux_effective)[0], flux_effective_func(flux_effective)[1], label='Effective flux', color='blue')
# ax2.plot(flux_theoretical[:, 0], flux_theoretical[:, 1], label='Effective flux', color='blue')

ax.set_title('Particle flux EP\(_1\)')
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Particle flux (particles/cm\(^2\)/s)')

# ax2.set_title('Increase in radius')
# ax2.set_xlabel('Time (ps)')
# ax2.set_ylabel('Radius (Å)')

df.plot(x='Time', y=['Atoms'], ax=ax2, color='blue', legend=False)
ax2.set_title('Number of particles EP\(_1\)')
ax2.set_xlabel('Time (ps)')
ax2.set_ylabel('Number of particles')
ax2.set_ylim(8.2784e6, 8.281e6)
# ax2.set_ylim(7.89e6, 8.3e6)

plt.tight_layout()
# plt.show()

plt.savefig('figs/flux_atoms_f25.pdf', dpi=300, bbox_inches='tight')
