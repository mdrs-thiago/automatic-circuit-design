import pygad 

def change_mutation(ga):
    if ga.generations_completed > 0.5 * ga.num_generations:
        new_mutation += ga.mutation_percent_genes * ga.multiply_mutation
        new_mutation = max(0.9, new_mutation)

class GeneticAlgorithmPyGAD:

    def __init__(self, ngenes: int, fitness_function, generations: int = 50, num_parents_mating: int = 2,
                popsize: int = 10, crossover_type: str = 'single_point', gene_space: dict = {},
                mutation_type: str = 'random', parent_selection: str = 'sss', gene_type: list = ['float'],
                mutation_rate: float = 0.01, crossover_rate: float = 0.7, multiply_mutation: float = None):

        self.ngenes = ngenes 
        self.fitness_function = fitness_function
        self.generations = generations
        self.popsize = popsize 
        self.gene_space = gene_space
        self.crossover_type = crossover_type
        self.mutation_type = mutation_type 
        self.gene_type = gene_type 
        self.parent_selection = parent_selection
        self.mutation_rate = mutation_rate
        self.crossover_probability = crossover_rate
        self.num_parents_mating = num_parents_mating
        
        on_generation = None
        if multiply_mutation is not None:
            on_generation = change_mutation

        self.ga_instance = pygad.GA(num_generations=self.generations,
                            num_parents_mating=self.num_parents_mating,
                            fitness_func=self.fitness_function,
                            sol_per_pop=self.popsize,
                            num_genes=self.ngenes,
                            #init_range_low=self.lb,
                            #init_range_high=self.ub,
                            gene_space = self.gene_space,
                            parent_selection_type=self.parent_selection,
                            #keep_parents=keep_parents,
                            crossover_type=self.crossover_type,
                            crossover_probability = self.crossover_probability,
                            mutation_type=self.mutation_type,
                            mutation_percent_genes=self.mutation_rate,
                            on_generation=on_generation,
                            gene_type = self.gene_type)

        self.ga_instance.multiply_mutation = multiply_mutation

    def run(self):
        self.ga_instance.run()


