import random


class GeneticAlgorithmUtils:
    system_random = random.SystemRandom()

    @classmethod
    def get_random_element_for_mutation(cls, valid_gene_data):
        return GeneticAlgorithmUtils.system_random.choice(valid_gene_data)

    @classmethod
    def get_random_number_from_list(cls, number_list):
        return GeneticAlgorithmUtils.system_random.choice(number_list)
