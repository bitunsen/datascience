import random
from Gene import Gene
from GeneticAlgorithmUtils import GeneticAlgorithmUtils


class Chromosome:
    _gene_array = []
    _point_array = []
    fitness_score = 0
    _valid_gene_data = []
    _mutate_ratio = 0.0

    def __init__(self, _data_array, valid_gene_data):
        self._gene_array = _data_array[:]
        self._valid_gene_data = valid_gene_data
        for i in range(len(_data_array)):
            _point = Gene(i, _data_array[i])
            self._point_array.append(_point)
        _mutate_ratio = (len(set(_data_array))/len(_data_array))

    def get_mutate_ratio(self):
        return self._mutate_ratio

    def print_chromosome(self, message):
        gene_size = len(self._gene_array)
        chromosome_str = "["
        for i in range(gene_size):
            if i < (gene_size - 1):
                chromosome_str = chromosome_str + str(self._gene_array[i]) + ", "
            if i == (gene_size - 1):
                chromosome_str = chromosome_str + str(self._gene_array[i])
        chromosome_str = chromosome_str + "]"
        print(message, "Chromosome : ", chromosome_str, " ::: Score : ", self.fitness_score)

    def get_fitness_score(self):
        return self.fitness_score

    def set_fitness_score(self, score):
        self.fitness_score = score

    def get_gene_array(self):
        return self._gene_array

    def get_all_elements_as_gene(self):
        return self._point_array

    def crossover(self, other_parent):
        index_crossover_array = [i for i in range(len(self.get_gene_array()))]
        index_crossover = GeneticAlgorithmUtils.get_random_number_from_list(index_crossover_array)
        if index_crossover == 0 or index_crossover == (len(self.get_gene_array()) - 1):
            index_crossover = int(len(self.get_gene_array()) / 2)
        # index_crossover = int(len(self.get_gene_array()) / 2)
        parent_1_data = self.get_gene_array()
        parent_2_data = other_parent.get_gene_array()
        child_data_1 = parent_1_data[0:index_crossover]
        child_data_2 = parent_2_data[0:index_crossover]
        child_data_1.extend(parent_2_data[index_crossover:])
        child_data_2.extend(parent_1_data[index_crossover:])

        child_1 = Chromosome(child_data_1, self._valid_gene_data)
        child_2 = Chromosome(child_data_2, self._valid_gene_data)

        # self.mutation(self._valid_gene_data)
        # other_parent.mutation(self._valid_gene_data)
        # print(child_2)
        return child_1, child_2

    def mutation(self):
        _random_sample_index_array = random.sample(range(len(self._gene_array)), 2)
        _index_1_of_mutation = _random_sample_index_array[0]
        _index_2_of_mutation = _random_sample_index_array[1]

        self._gene_array[_index_1_of_mutation], self._gene_array[_index_2_of_mutation] = self._gene_array[
                                                                                           _index_2_of_mutation], \
                                                                                       self._gene_array[
                                                                                           _index_1_of_mutation]
        point = Gene(_index_1_of_mutation, self._gene_array[_index_1_of_mutation])
        self._point_array[_index_1_of_mutation] = point
        point = Gene(_index_2_of_mutation, self._gene_array[_index_2_of_mutation])
        self._point_array[_index_2_of_mutation] = point
