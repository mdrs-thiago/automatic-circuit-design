from optimizer import GeneticAlgorithmPyGAD as GA
from ltspice_circuits import Circuit
from revalidate_netlist import Revalidator
import numpy as np 
import matplotlib.pyplot as plt 

def fitness_function_parametrizable(components:list, operations:list, analysis:str = 'get_data', operation:str = '.op V(2)', node:str='V(2)', signal_ref = 2.5, plot = False):
    def fitness_function(solution, solution_idx):
        r_values = [10, 11, 12, 13, 15, 16, 18, 20, 22, 24, 27, 30, 33, 36, 39, 43, 47, 51, 56, 62, 68, 75, 82, 91]
        common_components = ('L','C','Q','D', 'V')
        list_components = ['Optimized circuit \n']
        i = 0
        components_in_node = {}
        for component in components:
            if component[0][0] == 'R':
                comp = np.array(component)
                s = sum(comp == None)
                if comp[-1] is None:
                    #If the last value of resistor is empty, we use two chromosomes 
                    #to represent its value. 
                    comp_sol = solution[i:i+s+1].copy()
                    
                    comp_sol[-2] = r_values[comp_sol[-2]] * 10**comp_sol[-1]
                    comp_sol = comp_sol[:-1]
                    __component = write_component(comp, comp_sol)
                    i += s + 1
                elif s == 0:
                    __component = ' '.join(component)
                
                else:
                    comp_sol = solution[i:i+s]
                    __component = write_component(comp, comp_sol)
                    i += s
                
            elif component[0].startswith(common_components):
                comp = np.array(component)
                s = sum(comp == None)
                if s == 0:
                    __component = ' '.join(component)
                else:
                    comp_sol = solution[i:i+s]
                    __component = write_component(comp, comp_sol)
                    i += s
                __comp_sliced = __component.split()
                if __comp_sliced[0][0] == 'Q':
                    for terminal in __comp_sliced[1:4]:
                        if terminal not in components_in_node.keys():   
                            components_in_node[terminal] = 1
                        else:
                            components_in_node[terminal] += 1 
                else:
                    for terminal in __comp_sliced[1:3]:
                        if terminal not in components_in_node.keys():   
                            components_in_node[terminal] = 1
                        else:
                            components_in_node[terminal] += 1
            else:
                #print(component)
                __component = ' '.join(component)

            #Adicionando capacitores com uF
            if component[0][0] == 'C':
                __component += 'u'

            list_components.append(__component)
        
        

        #Revalidação do circuito
        revalidator = Revalidator(list_components)
        revalidator.revalidate(components_in_node) 
        list_components = revalidator.list_components

        list_components.extend(operations)

        netcir = Circuit(analysis) 
        netcir.generate_cir(list_components) 
        netcir.run_cir()

        #Transformar essa parte em uma nova função para ajustar aos mais 'comuns'.
        type_op = operation.split()[0]
        try:
            if 'op' in type_op:
                res = netcir.read_output(node)
                return op_diff(res[0], signal_ref)
            elif 'tran' in type_op or 'dc' in type_op:
                res, t = netcir.read_output(node)
                #signal_ref = expected_signal(t)
                if plot:
                    plt.plot(res[1:])
                    plt.plot(signal_ref)
                    plt.show()
                
                diff = abs(res[1:] - signal_ref) 
                #return -np.sum(diff) - 10*np.sum(diff[40:80])
                return vector_diff(res[1:], signal_ref)
                
        except Exception as e:
            print(e)
            return -10e9
        
    return fitness_function

def vector_diff(res, ans):
    return -np.mean(abs(res - ans))


def op_diff(res: float, ans: float):
    '''
    Diferença absoluta em um ponto de operação 
    '''
    return -abs(res - ans)


def write_component(component, sol):
    c = np.array(component) 
    c[c == None] = sol
    netcomp = c.astype('str').tolist()
    netcomp = ' '.join(netcomp)
    
    return netcomp

def expected_signal(t):
    return 4.5*np.sin(2*np.pi*300*t + np.pi) + 4.5

def test_resistive():
    components = [['V1','1','0','5'],['R1','1','2', None],['R2','2','0', None]]
    operations = ['.op V(2)', '.end']
    fitness_function = fitness_function_parametrizable(components, analysis='get_data',operation='.op V(2)', node='V(2)')

    ga = GA(ngenes = 2, fitness_function = fitness_function, 
            gene_space=[{'low':100,'high':1000},{'low':1000,'high':10000}])
    ga.run()

def test_amplifier(): 
    components = [['Q1', 'Out', None, None, '0', 'NPN'], ['R1', None, 'Out', None], ['R2', None, None, None], 
                  ['V1', None, '0', '9'], ['V2', '5', '0', 'SINE(0 0.05 300)'], ['R3', None, None, None],
                  ['R4', None, None, None], ['C1', None, None, None], ['C2', None, None, None]]

    operations = ['.model NPN NPN', '.model PNP PNP', 
                  '.lib C:\\Users\\thiag\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.bjt',
                  '.tran 1m 500m', '.backanno', '.end']
    t = {'low': 0, 'high':8}
    r = {'low': 0, 'high': 23}
    p = {'low': 1, 'high': 7}
    c = {'low':0.01, 'high':1000}
    gene_space = [t, t, t, r, t, t, r, t, t, t, r, t, t, r, t, t, c, t, t, c]
    gene_type = [int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, [float,2], int, int, [float,2]]

    fitness_function = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                                       operation='.tran 1m 500m', node='V(Out)')

    ga = GA(ngenes = len(gene_type), fitness_function = fitness_function, gene_type = gene_type, gene_space=gene_space, popsize=5, generations=5)
    best_solution = ga.run()
    print('Running best circuit')
    fitness_function(best_solution, None)


def test_sigmoid(): 
    components = [['R1', 'N001', 'N002', None], ['R2', 'N001', 'Out', None], 
                  ['R3', 'N006', 'N008', None], ['R4', 'N005', 'N008', None],
                  ['R5', 'N008', 'N007', None], ['Q1', 'N002', 'N003', 'N006', '0', '2N3904'],
                  ['Q2', 'Out', 'N004', 'N005', '0', '2N3904'], ['V1', 'N003', '0', '0'],
                  ['V2', 'N004', '0', None], ['V3', 'N001', '0', None], ['V4', 'N007', '0', None]] 

    operations = ['.model NPN NPN', '.model PNP PNP', 
                  '.lib C:\\Users\\thiag\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.bjt',
                  '.dc V1 -3 3 0.05', '.backanno', '.end']
    
    t = {'low': 0, 'high':10}
    r = {'low': 100, 'high': 10000}
    c = {'low':0.01, 'high':1000}
    v = {'low':-15, 'high':15}
    gene_space = [r, r, r, r, r, v, v, v]
    gene_type = [int]*len(gene_space)

    v_in = np.arange(-3, 3, 0.05)
    v_out = 1/(1 + np.exp(-6*v_in))

    fitness_function = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                                       operation='.dc -3 3 0.05', node='V(Out)', signal_ref = v_out)

    ga = GA(ngenes = len(gene_type), parent_selection='sss', fitness_function = fitness_function, gene_type = gene_type, gene_space=gene_space, popsize=50, generations=50, mutation_rate=0.2)
    best_solution = ga.run()
    print('Running best circuit')
    best_ = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                            operation='.dc -3 3 0.05', node='V(Out)', signal_ref = v_out, 
                                            plot=True)

    best_(best_solution, None)

def test_soft_clipping(): 
    components = [['D1', 'N002', 'N006', 'D'], ['D2', 'N002', 'N005', 'D'], ['D3', 'N004', 'N002', 'D'],
                  ['D4', 'N003', 'N002', 'D'], ['R1', 'out', 'N003', None], ['R2', 'out', 'N004', None],
                  ['R3', 'out', 'N005', None], ['R4', 'out', 'N006', None], 
                  ['XU1', 'N009', 'N002', 'N007', 'N010', 'out', 'AD549'], ['R5', 'N009', '0', None],
                  ['V1', 'N010', '0', None], ['V2', 'N007', '0', None], ['R6', 'N002', '0', None],
                  ['R7','N002','out',None], 
                  #['R7', 'N006', 'N001', None], ['R8', 'N005', 'N001', None], ['V4', 'N001', '0', None],
                  ['V3', 'N008', '0', '0'], ['R9', 'N009', 'N008', None]]

    operations = ['.model D D', '.lib C:\\Users\\thiag\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio',
                  '.dc V3 -5 5 0.1', '.lib ADI1.lib', '.backanno', '.end']
    
    t = {'low': 0, 'high':7}
    r = {'low': 100, 'high': 10000}
    c = {'low':0.01, 'high':1000}
    v = {'low':-15, 'high':15}
    gene_space = [r, r, r, r, r, v, v, r, r, r]
    gene_type = [int]*len(gene_space)

    v_in = np.arange(-5, 5, 0.1)
    v_out = 1/(1 + np.exp(-6*v_in))

    fitness_function = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                                       operation='.dc -5 5 0.1', node='V(out)', signal_ref = v_out)

    ga = GA(ngenes = len(gene_type), parent_selection='sss', fitness_function = fitness_function, gene_type = gene_type, gene_space=gene_space, popsize=50, generations=50, mutation_rate=0.2)
    best_solution = ga.run()
    print('Running best circuit')
    best_ = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                            operation='.dc -5 5 0.1', node='V(out)', signal_ref = v_out, 
                                            plot=True)

    best_(best_solution, None)

'''def test_s_mf_function(): 
    components = [['R2', '1', '5', None], ['R3', '5', '0', None], ['R4', '5', '6', None], 
                  ['XU1', '4', '5', '6', 'N010', 'N007', 'LT1803'], ['R5', None, None, None],
                  ['V1', 'N010', '0', None], ['V2', 'N007', '0', None],
                  #['R7','N002','out',None], 
                  #['R7', 'N006', 'N001', None], ['R8', 'N005', 'N001', None], ['V4', 'N001', '0', None],
                  ['V3', '4', '0', '0'], ['V4', '1', '0', None]]

    
    
# components = [['D1', None, None, 'D'], ['D2', None, None, 'D'], ['R1', None, None, None], 
#                   ['R2', None, None, None], ['R3', None, None, None], ['R4', None, None, None], 
#                   ['XU1', '4', '5', 'N007', 'N010', '6', 'AD549'], ['R5', None, None, None],
#                   ['V1', 'N010', '0', None], ['V2', 'N007', '0', None], ['R6', None, None, None],
#                   #['R7','N002','out',None], 
#                   #['R7', 'N006', 'N001', None], ['R8', 'N005', 'N001', None], ['V4', 'N001', '0', None],
#                   ['V3', None, '0', '0'], ['V4', None, '0', None]]


    operations = ['.model D D', '.lib C:\\Users\\thiag\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio',
                  '.dc V3 -2 1 0.1', '.lib LTC2.LIB', '.backanno', '.end']
    
    t = {'low': 0, 'high':7}
    r = {'low': 0, 'high': 23}
    p = {'low': 1, 'high': 7}
    c = {'low':0.01, 'high':1000}
    v = {'low':-15, 'high':15}
    #gene_space = [t, t, t, t, t, t, r, t, t, r, t, t, r, t, t, r, t, t, r, v, v, t, t, r, t, t, v]
    gene_space = [r, p, r, p, r, p, t, t, r, p, v, v, v]

    gene_type = [int]*len(gene_space)
    print(len(gene_space))
    v_1 = np.zeros(10)
    v_2 = np.arange(-0.9, 0.1, 0.1) + 1
    v_3 = np.ones(10)
    v_out = np.hstack((v_1, v_2, v_3))
    fitness_function = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                                       operation='.dc -2 1 0.1', node='V(6)', signal_ref = v_out)

    ga = GA(ngenes = len(gene_type), parent_selection='sss', fitness_function = fitness_function, gene_type = gene_type, gene_space=gene_space, popsize=50, generations=50, mutation_rate=0.25)
    best_solution = ga.run()
    print('Running best circuit')
    best_ = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                            operation='.dc -2 1 0.1', node='V(6)', signal_ref = v_out, 
                                            plot=True)

    best_(best_solution, None)'''


def test_s_mf_function(): 
    components = [['R1', '1', '3', None], ['R3', '3', '0', None], ['R4', '3', '6', None], 
                  ['XU1', '4', '3', '6', 'N010', 'N007', 'LT1803'], ['R5', None, None, None],
                  ['V1', 'N010', '0', None], ['V2', 'N007', '0', None],
                  #['R7','N002','out',None], 
                  #['R7', 'N006', 'N001', None], ['R8', 'N005', 'N001', None], ['V4', 'N001', '0', None],
                  ['V3', '4', '0', '0'], ['V4', '1', '0', None]]

# components = [['D1', None, None, 'D'], ['D2', None, None, 'D'], ['R1', None, None, None], 
#                   ['R2', None, None, None], ['R3', None, None, None], ['R4', None, None, None], 
#                   ['XU1', '4', '5', 'N007', 'N010', '6', 'AD549'], ['R5', None, None, None],
#                   ['V1', 'N010', '0', None], ['V2', 'N007', '0', None], ['R6', None, None, None],
#                   #['R7','N002','out',None], 
#                   #['R7', 'N006', 'N001', None], ['R8', 'N005', 'N001', None], ['V4', 'N001', '0', None],
#                   ['V3', None, '0', '0'], ['V4', None, '0', None]]


    operations = ['.model D D', '.lib C:\\Users\\thiag\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.dio',
                  '.dc V3 -2 1 0.1', '.lib LTC2.LIB', '.backanno', '.end']
    
    t = {'low': 0, 'high':7}
    r = {'low': 0, 'high': 23}
    p = {'low': 1, 'high': 7}
    c = {'low':0.01, 'high':1000}
    v = {'low':-15, 'high':15}
    #gene_space = [t, t, t, t, t, t, r, t, t, r, t, t, r, t, t, r, t, t, r, v, v, t, t, r, t, t, v]
    gene_space = [r, p, r, p, r, p, t, t, r, p, v, v, v]

    gene_type = [int]*len(gene_space)
    print(len(gene_space))
    v_1 = np.zeros(10)
    v_2 = np.arange(-0.9, 0.1, 0.1) + 1
    v_3 = np.ones(10)
    v_out = np.hstack((v_1, v_2, v_3))
    fitness_function = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                                       operation='.dc -2 1 0.1', node='V(6)', signal_ref = v_out)

    ga = GA(ngenes = len(gene_type), parent_selection='sss', fitness_function = fitness_function, gene_type = gene_type, gene_space=gene_space, popsize=50, generations=150, mutation_rate=0.25)
    best_solution = ga.run()
    print('Running best circuit')
    best_ = fitness_function_parametrizable(components, operations, analysis='ac_get_data', 
                                            operation='.dc -2 1 0.1', node='V(6)', signal_ref = v_out, 
                                            plot=True)

    best_(best_solution, None)







if __name__ == '__main__':
    test_s_mf_function()