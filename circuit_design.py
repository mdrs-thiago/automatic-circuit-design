from optimizer import GeneticAlgorithmPyGAD as GA
from ltspice_circuits import Circuit
import numpy as np 

def fitness_function_parametrizable(components:list, operation:str = '.op V(2)', node:str='V(2)', signal_ref = 2.5):
    def fitness_function(solution, solution_idx):
        common_components = ('R','L','C')
        list_components = ['Optimized circuit \n']
        i = 0
        for component in components:
            if component[0].startswith(common_components):
                s = 4 - len(component)
                comp_sol = solution[i:i+s]
                __component = write_component(component, comp_sol)
                i += s
            else:
                __component = ' '.join(component)

            list_components.append(__component)

        list_components.append(f'{operation}\n')
        list_components.append('.end')
        netcir = Circuit('get_data') 
        netcir.generate_cir(list_components) 
        netcir.run_cir()

        #Transformar essa parte em uma nova função para ajustar aos mais 'comuns'.
        
        res = netcir.read_output(node)
        type_op = operation.split()[0]
        if type_op.contains('op'):
            return op_diff(res[0], signal_ref)
        elif type_op.contais('tran'):
            return vector_diff(res, signal_ref)
        
    return fitness_function

def vector_diff(res, ans):
    return -np.mean(abs(res - ans))


def op_diff(res: float, ans: float):
    '''
    Diferença absoluta em um ponto de operação 
    '''
    return -abs(res - ans)


def write_component(component, sol):
    s = len(component)
    if s == 1:
        netcomp = f'{component[0]} {sol[0]} {sol[1]} {sol[2]}'
    elif s == 2:
        netcomp = f'{component[0]} {sol[0]} {sol[1]} {component[1]}'
    elif s == 3:
        netcomp = f'{component[0]} {component[1]} {component[2]} {sol[0]}'
    else:
        netcomp = f'{component[0]} {component[1]} {component[2]} {component[3]}'
    
    return netcomp


if __name__ == '__main__':
    components = [('V1','1','0','5'),('R1','1','2'),('R2','2','0')]

    fitness_function = fitness_function_parametrizable(components, operation='.op V(2)', node='V(2)')

    ga = GA(ngenes = 2, fitness_function = fitness_function, 
            gene_space=[{'low':100,'high':1000},{'low':1000,'high':10000}])
    ga.run()