import numpy as np
import random
from Chromosome import Chromosome
from Gene import Gene


def get_all_attacking_points_from_one_point(x_cor, y_cor, num_of_queens):
    _invalid_point_arrays = []
    for i in range(num_of_queens):
        if i == x_cor:
            continue
        offset = (x_cor - i)

        if 0 <= (y_cor + offset) < num_of_queens:
            value = str(i) + "-" + str(y_cor + offset)
            _invalid_point_arrays.append(Gene(i, (y_cor + offset)))

        if 0 <= (y_cor - offset) < num_of_queens:
            _invalid_point_arrays.append(Gene(i, (y_cor - offset)))

    for i in range(num_of_queens):
        if i == x_cor:
            continue
        _invalid_point_arrays.append(Gene(i, y_cor))
        # print("Point :: (",x_cor,", ",i,")")
    return _invalid_point_arrays


class GAExecution:
    _attacking_point_dictionary = {}
    _size_of_population = 0
    _length_of_chromosome = 0
    _valid_gene_data = []
    _population_pool = []
    _mutation_rate = 0.01
    _max_generation = 2000
    _system_random = random.SystemRandom()
    _current_generation = 0

    def __init__(self, size_of_population, num_of_queens, mutation_rate):
        self._size_of_population = size_of_population

        self._mutation_rate = mutation_rate

        self._length_of_chromosome = num_of_queens
        self._valid_gene_data = [i for i in range(num_of_queens)]
        self.__populate_attacking_point_dictionary()

    def get_current_generation(self):
        return self._current_generation

    def set_current_generation(self, cur_generation):
        self._current_generation = cur_generation

    def __populate_attacking_point_dictionary(self):
        num_of_queens = self._length_of_chromosome
        for row in range(num_of_queens):
            for column in range(num_of_queens):
                key = str(row) + "-" + str(column)
                self._attacking_point_dictionary[key] = get_all_attacking_points_from_one_point(row, column,
                                                                                                num_of_queens)

    def get_all_attacking_points_from_dictionary(self, gene_element):
        x_cor, y_cor = gene_element.x_cor, gene_element.y_cor
        key = str(x_cor) + "-" + str(y_cor)
        return self._attacking_point_dictionary[key]

    def generate_initial_population(self):

        self._current_generation = 1
        for count_population in range(self._size_of_population):
            # print("Valid Data Array :: ", self._valid_gene_data)
            gene_array = self._valid_gene_data[:]
            random.shuffle(gene_array)
            # print("Gene Array :: ", gene_array)
            _initial_chromosome = Chromosome(gene_array, self._valid_gene_data)
            _chromosome_score = self.calculate_fitness_score(_initial_chromosome)
            self._population_pool.append(_initial_chromosome)

        # Sorting (in ascending order) all the chromosomes present in the population pool
        self._population_pool = sorted(self._population_pool, key=lambda element: element.fitness_score,
                                       reverse=True)
        # for _element in self._population_pool:
        #    _element.print_chromosome("Initial Pool :: ")

        # print("-------------- END OF POPULATION GENERATION AND SORTING BASED ON SCORE --------------------------")
        return self._population_pool

    def get_population_pool_size(self):
        return len(self._population_pool)

    def calculate_fitness_score(self, input_chromosome):
        """

        :type input_chromosome: Chromosome
        """
        _data = input_chromosome.get_gene_array()[:]
        _non_attacking_queen = [0, 0, 0, 0, 0, 0, 0, 0]
        for index in range(len(_data)):
            current_point = Gene(index, _data[index])
            # current_point.print("Current Queen :: ")
            attacking_queen = 0
            current_attacking_points = self.get_all_attacking_points_from_dictionary(current_point)
            for i in range(len(_data)):
                if i == index:
                    continue
                other_queen = Gene(i, _data[i])
                for attacking_point in current_attacking_points:
                    if other_queen.equals(attacking_point):
                        attacking_queen = attacking_queen + 1
                        # other_queen.print("Attacking Queen :: ")

            _non_attacking_queen[index] = ((len(_data) - 1) - attacking_queen)
            # print(non_attacking_queen)
            input_chromosome.fitness_score = np.sum(_non_attacking_queen)
        return input_chromosome.fitness_score

    # Take half of the high scored parents for producing next generation kids
    def run_parent_selection_process(self):
        _max_score = 0
        _score_array = []
        _selected_parent_pool = []
        _selected_parent_index_set = set()
        initial_pool_size = len(self._population_pool)
        _selected_parent_pool_size = int(self._size_of_population / 2)

        if _selected_parent_pool_size % 2 != 0:
            _selected_parent_pool_size += 1

        _selected_parent_pool = self._population_pool[0:_selected_parent_pool_size]
        return _selected_parent_pool, _max_score

    # Select parents using tournament process
    # Using Tournament process, Only half of the parents will be selected for reproduction
    def run_parent_selection_process_tournament(self):
        _max_score = 0
        _score_array = []
        _selected_parent_pool = []
        _selected_parent_index_set = set()
        initial_pool_size = len(self._population_pool)
        _selected_parent_pool_size = int(self._size_of_population / 2)

        # print("-------------- START OF TOURNAMENT METHOD OF PARENT SELECTION --------------------------")
        # print("[GEN - ", self._current_generation, "] :: Number of parents to be selected :: ",
        #      _selected_parent_pool_size)

        # Calculate Fitness Score (Number of non-attacking Queens) for each chromosome
        for _chromosome in self._population_pool:
            # _chromosome_score = self.calculate_fitness_score(_chromosome)
            _score_array.append(_chromosome.get_fitness_score())

        if len(self._population_pool) > self._size_of_population:
            self._population_pool = self._population_pool[:self._size_of_population]

        _max_score = np.max(_score_array)

        _count = 0
        while len(_selected_parent_index_set) < _selected_parent_pool_size:
            _random_sample_index = random.sample(range(len(self._population_pool)), 4)
            _count += 1
            if _count % 10 == 0:
                print("Infinite Loop ....")
            _max_selected_score = -1
            _max_selected_index = -1
            for _index in _random_sample_index:
                if self._population_pool[_index].get_fitness_score() > _max_selected_score:
                    _max_selected_index = _index
            if _max_selected_index in _selected_parent_index_set:
                continue
            else:
                _selected_parent_pool.append(self._population_pool[_max_selected_index])
                _selected_parent_index_set.add(_max_selected_index)

        # for _chromosome in _selected_parent_pool:
        #    _chromosome.print_chromosome("[GEN - " + str(self._current_generation) + "] :: Selected Parent :: ")

        # selected_parent_pool = self._population_pool[:]

        return _selected_parent_pool, _max_score

    def execute_crossover_mutation(self, _selected_parent_pool):
        _num_of_selected_parents = len(_selected_parent_pool)
        _mating_pool = []
        # print("-------------- START OF CROSSOVER AND MUTATION OF PARENT SELECTION --------------------------")
        for _index in range(_num_of_selected_parents):
            _index_1 = _index % _num_of_selected_parents
            _index_2 = (_index + 1) % _num_of_selected_parents
            _mating_pool.append(_selected_parent_pool[_index_1])
            _mating_pool.append(_selected_parent_pool[_index_2])

        # for _chromosome in _mating_pool:
        #    _chromosome.print_chromosome("[GEN - " + str(self._current_generation) + "] :: Matting Parent :: ")

        self._population_pool = []
        for _index in range(int(len(_mating_pool) / 2)):
            _parent_1 = _mating_pool[(2 * _index)]
            _parent_2 = _mating_pool[(2 * _index) + 1]
            _child_1, _child_2 = _parent_1.crossover(_parent_2)

            if _child_1.get_mutate_ratio() < 0.6:
                _child_1.mutation()
            if _child_2.get_mutate_ratio() < 0.6:
                _child_2.mutation()

            self.calculate_fitness_score(_child_1)
            self.calculate_fitness_score(_child_2)
            self._population_pool.append(_child_1)
            self._population_pool.append(_child_2)

        # print("-------------- AFTER CROSSOVER AND MUTATION AND BEFORE SORTING OF "
        #      "PARENT SELECTION --------------------------")

        # for _chromosome in self._population_pool:
        #    _chromosome.print_chromosome("[GEN - " + str(self._current_generation)
        #    + "] :: Before Sorting New Child :: ")

        # Sorting (in ascending order) all the chromosomes present in the population pool
        self._population_pool = sorted(self._population_pool, key=lambda element: element.fitness_score,
                                       reverse=True)
        # print("-------------- AFTER CROSSOVER AND MUTATION AND SORTING OF PARENT SELECTION --------------------------")

        # for _chromosome in self._population_pool:
        #    _chromosome.print_chromosome("[GEN - " + str(self._current_generation) + "] :: New Child :: ")

        return self._population_pool


if __name__ == "__main__":

    valid_data = [i for i in range(8)]
    engine = GAExecution(100, 8, 0.01)
    current_generation = 1
    _max_generation = 2000
    target_chromosome_array = [Chromosome([0, 0, 0, 0, 0, 0, 0, 0], valid_data)]
    _initial_population_pool = engine.generate_initial_population()
    _max_scored_chromosome = _initial_population_pool[0]

    if _max_scored_chromosome.get_fitness_score() == 56:
        target_chromosome_array[0] = _max_scored_chromosome
        _max_scored_chromosome.print_chromosome(
            "Generation #" + str(current_generation) + " Target Chromosome :: ")
    else:

        while True:
            selected_parent_pool, max_score = engine.run_parent_selection_process()
            _initial_population_pool = engine.execute_crossover_mutation(selected_parent_pool)
            _max_scored_chromosome = _initial_population_pool[0]
            current_generation = current_generation + 1
            engine.set_current_generation(current_generation)

            if _max_scored_chromosome.get_fitness_score() == 56:
                target_chromosome_array[0] = _max_scored_chromosome
                _max_scored_chromosome.print_chromosome(
                    "Generation #" + str(current_generation) + " Target Chromosome :: ")
                break

            if current_generation > _max_generation:
                print("Number of Maximum Generation limit crossed. Terminating the engine....")
                break

            highest_score_candidate = _initial_population_pool[0]
            second_highest_score_candidate = _initial_population_pool[1]
            if current_generation % 10 == 0:
                highest_score_candidate.print_chromosome("Generation #" + str(current_generation) + " ::: Highest ")
                second_highest_score_candidate.print_chromosome(
                    "Generation #" + str(current_generation) + " ::: Second Highest ")
