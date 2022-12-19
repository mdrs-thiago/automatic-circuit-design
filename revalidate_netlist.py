import numpy as np 
from copy import copy

class Revalidator:
    def __init__(self, list_components:list):
        self.list_components = list_components

    def check_floating_nodes(self,components_in_node:list):
        for k,v in components_in_node.items():
            if v == 1:
                self.list_components.append(f'R99{k} {k} 0 10')

    def check_short_circuit_transistor(self):
        for v in self.list_components:
            if v.startswith('Q'):
                __v = np.array(v.split())
                terminals, counts = np.unique(__v[1:4], return_counts=True)
                sct = {t:c for t,c in zip(terminals, counts) if c > 1}

                if len(sct) > 0:
                    new_list = copy(self.list_components)   
                    for t,c in sct.items():
                        
                        while c > 1:
                            new_terminal = f'N{t}{str(c)}'
                            __v[np.argwhere(__v==t)[0]] = new_terminal
                            new_list = [' '.join(__v) if k == v else k for k in new_list]
                            new_list.append(f'R{__v[0]}{c} {new_terminal} 0 1000')
                            c -= 1

                    self.list_components = new_list

    def check_parallel_sources(self):
        v_sources = np.array([c for c in self.list_components if c.startswith('V')])
        v_terminals = np.array([' '.join(c.split()[1:3]) for c in self.list_components if c.startswith('V')])
        
        if len(v_terminals) > 1:
            new_list = copy(self.list_components)
            mask = np.array([True if np.sum(t == v_terminals) > 1 else False for t in v_terminals])
            for i,v in enumerate(v_sources[mask]):
                v_components = v.split()
                new_resistor = f'R9{i} N10{i}{v_components[1]} {v_components[1]} 1000'
                v_components[1] = f'N10{i}{v_components[1]}'
                new_list.append(new_resistor)
                new_list = [' '.join(v_components) if k == v else k for k in new_list]
            
        self.list_components = new_list

    def revalidate(self, components_in_node):
        self.check_short_circuit_transistor()
        self.check_floating_nodes(components_in_node)
        self.check_parallel_sources()
                
                
        

if __name__ == '__main__':
    list_components =  ['Q1 Out 5 5 0 NPN', 'R1 6 Out 56471', 'R2 0 5 62740',
                        'V1 5 0 9', 'V2 5 0 SINE(0 0.05 300)', 'R3 3 0 91614',
                        'R4 7 0 55723', 'C1 6 3 6.68u', 'C2 3 0 353.58u']
    
    operations = ['.model NPN NPN', '.model PNP PNP',
                  '.lib C:\\Users\\thiag\\OneDrive\\Documents\\LTspiceXVII\\lib\\cmp\\standard.bjt',
                  '.tran 1m 500m', '.backanno', '.end']

    components_in_node = {}
    for component in list_components:
        __comp_sliced = component.split()
        if __comp_sliced[0].startswith('Q'):
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

    revalidate = Revalidator(list_components)

    revalidate.revalidate(components_in_node)
    print(revalidate.list_components)