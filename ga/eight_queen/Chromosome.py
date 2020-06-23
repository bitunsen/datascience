import random
from Gene import Gene
from GeneticAlgorithmUtils import GeneticAlgorithmUtils

"""
Chromosome is a class containing an array of Gene
Mutation rate is calculated based on the uniqueness of data elements in the data array
It is calculated by dividing the number of unique elements by the length of data elements
 
"""


class Chromosome:
    _data_array = []
    _gene_array = []
    fitness_score = 0
    _mutate_ratio = 0.0

    def __init__(self, _input_data_array):
        if not isinstance(_input_data_array, list):
            raise TypeError("Invalid argument passed - only array of positive integers are allowed.")

        if len(_input_data_array) < 2:
            raise ValueError("Invalid argument passed - number of elements must be more than 2 in the input array.")

        for _element in _input_data_array:
            if not isinstance(_element, int):
                raise TypeError("Invalid argument passed - only array of positive integers are allowed.")
            if _element < 0:
                raise TypeError("Invalid argument passed - only array of positive integers are allowed.")

        self._data_array = _input_data_array[:]
        for i in range(len(_input_data_array)):
            _point = Gene(i, _input_data_array[i])
            self._gene_array.append(_point)
        _mutate_ratio = (len(set(_input_data_array)) / len(_input_data_array))

    def get_mutate_ratio(self):
        return self._mutate_ratio

    def print_chromosome(self, message):
        gene_size = len(self._data_array)
        chromosome_str = "["
        for i in range(gene_size):
            if i < (gene_size - 1):
                chromosome_str = chromosome_str + str(self._data_array[i]) + ", "
            if i == (gene_size - 1):
                chromosome_str = chromosome_str + str(self._data_array[i])
        chromosome_str = chromosome_str + "]"
        print(message, "Chromosome : ", chromosome_str, " ::: Score : ", self.fitness_score)

    def get_fitness_score(self):
        return self.fitness_score

    def set_fitness_score(self, score):
        self.fitness_score = score

    def get_gene_array(self):
        return self._data_array

    def get_all_elements_as_gene(self):
        return self._gene_array

    """
    Rather than having fixed point crossover for all the parent chromosome, select 
    the crossover point randomly. This introduces more randomness in the chromosome during crossover 
    """
    def crossover(self, other_parent):
        if not isinstance(other_parent, Chromosome):
            raise TypeError("Invalid argument passed - only object of type Chromosome is allowed.")

        if len(other_parent.get_gene_array()) != len(self.get_gene_array()):
            raise TypeError(
                "Invalid argument passed - chromosome passed as argument does not have same number of Gene.")

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

        child_1 = Chromosome(child_data_1)
        child_2 = Chromosome(child_data_2)
        return child_1, child_2

    def mutation(self):
        if len(self._data_array) == 2:
            _index_1_of_mutation = 1
            _index_2_of_mutation = 0
        else:
            _random_sample_index_array = random.sample(range(len(self._data_array)), 2)
            _index_1_of_mutation = _random_sample_index_array[0]
            _index_2_of_mutation = _random_sample_index_array[1]

        self._data_array[_index_1_of_mutation], self._data_array[_index_2_of_mutation] = self._data_array[
                                                                                             _index_2_of_mutation], \
                                                                                         self._data_array[
                                                                                             _index_1_of_mutation]
        point = Gene(_index_1_of_mutation, self._data_array[_index_1_of_mutation])
        self._gene_array[_index_1_of_mutation] = point
        point = Gene(_index_2_of_mutation, self._data_array[_index_2_of_mutation])
        self._gene_array[_index_2_of_mutation] = point
