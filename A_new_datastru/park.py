def phase_transition(N_aminoacids, seq):
    A, J = initialize(N_aminoacids, seq)

    N_steps = 2000000
    ave_int = 1800000
    temp_int = 30
    temp_final = 20
    dtemp = 0.5
    temp_step = int((temp_int - temp_final) / dtemp) + 1
    temp_array = np.linspace(temp_int, temp_final, temp_step, endpoint=True)

    for t in temp_array:
        fname = 'UniJ5_equilibrium_{}_T{:.1f}.txt'.format(N_aminoacids, t)
        with open(fname, 'w') as file:
            for n in range(10):
                step, FreeEnergy, L, R = ProteinFolding_jit(N_aminoacids, J, A, t, N_steps)
                for i in range(len(step)):
                    if step[i] < ave_int:
                        continue
                    file.write('{} {} {}\n'.format(step[i], FreeEnergy[i], L[i]))

def equilibrium_mean_SD(fname):
    a = np.loadtxt(fname)
    E = a[:, 1]
    L = a[:, 2]
    std_E = np.std(E)
    std_L = np.std(L)
    se_E = std_E / np.sqrt(len(E))
    se_L = std_L / np.sqrt(len(L))
    return np.mean(E), np.mean(L), std_E, se_E, std_L, se_L

def phase_data(N_aminoacids=15):
    temp_int = 30
    temp_final = 0.5
    dtemp = 0.5
    temp_step = int((temp_int - temp_final) / dtemp) + 1
    temp_array = np.linspace(temp_int, temp_final, temp_step, endpoint=True)
    
    with open('UniJ5_phase_{}.txt'.format(N_aminoacids), 'w') as file:
        for t in temp_array:
            fname = "UniJ5_equilibrium_{}_T{:.1f}.txt".format(N_aminoacids, t)
            mean_E, mean_L, std_E, se_E, std_L, se_L = equilibrium_mean_SD(fname)
            file.write("{} {} {} {} {} {} {}\n".format(
                t, mean_E, se_E, std_E, mean_L, se_L, std_L))

