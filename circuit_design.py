from optimizer import GeneticAlgorithmPyGAD as GA
from ltspice_circuits import Circuit

def fitness_function_parametrizable(components):
    def fitness_function(solution, solution_idx):
        common_components = ('R','L','C')
        list_components = []
        i = 0
        for component in components:
            if component[0].startswith(common_components):
                s = 4 - len(component)
                comp_sol = solution[i:i+s]
                __component = write_component(component, comp_sol)
                i += s + 1
            else:
                __component = ' '.join(component)

            list_components.append(__component)

        netcir = Circuit('get_data') 
        netcir.generate_cir(list_components) 
        netcir.run_cir()

        #Transformar essa parte em uma nova função para ajustar aos mais 'comuns'.
        res = netcir.read_output('V(2)')
        return op_diff(res, 2.5)
    return fitness_function


def op_diff(res: float, ans: float):
    '''
    Diferença absoluta em um ponto de operação 
    '''
    return abs(res - ans)



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
    GA()
