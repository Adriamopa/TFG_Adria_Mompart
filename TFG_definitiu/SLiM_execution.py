from subprocess import Popen
import time


def repeated_iHS_SLIM(variable, path="/SLiM_recipe/calcul_iHS", output='/TFG_definitiu/vcfs_tfg/iHS', Q=1,
                         N=10000,sweep=0.7, L='1e6', S=0.01, rho_theta=500, MR=1, RR=1, num_rep=5):
    """
    This function receives at least a variable and runs the SLiM recipe in path with the variable values predefined 
    as parameters of the function plus the variable that has been chosen. The SLiM recipe is run num_rep number of times. 
    The output VCF files are saved in path. Possible variables that can be changed are: N (size of population), 
    sweep (sweep mutation frequency), L (chromosome length), S (selection coefficient), MR (mutation rate), RR (recombination rate).
    """
    dict = {'Q': Q, 'N': N, 'sweep': sweep, 'L': L, 'S': S, 'rho_theta': rho_theta, 'MR': MR, 'RR': RR, '': 'null'}
    print('SLiM iHS:')
    a = time.perf_counter_ns()
    reps = []
    for element in range(1,int(num_rep)+1):
        reps.append('slim {0}{1}{2}{3}{4}{5}{6}{7}{8}{9}'.format('-d Q='+str(Q), ' -d N='+str(N), ' -d sweep_freq='+str(sweep+0.1),
                                                           ' -d L='+str(L),' -d S='+str(S),' -d rho_theta='+str(rho_theta),
                                                           ' -d MR='+str(MR),' -d RR='+str(RR),' -d output='+"'"+output+'_'
                                                           +str(variable)+'_'+str(dict[variable])+'_'+str(element)+'.vcf'+"'", ' '+path))
    print(reps)
    processes = [Popen(cmd, shell=True) for cmd in reps]
    for p in processes: p.wait()
    b = time.perf_counter_ns()
    print(f"Total SLiM iHS execution time: {(b- a) * 1e-9}")




