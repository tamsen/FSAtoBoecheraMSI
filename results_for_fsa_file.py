import accuracy

class FSA_File_Results:
    def __init__(self, MSI_loci_results_by_loci):
        self.MSI_loci_results_by_loci = MSI_loci_results_by_loci


class MSI_loci_results:

    def __init__(self, alleles_called,
                 original_file_name, trace_data, ladder_data, mapping_data):
        self.alleles_called = alleles_called
        self.original_file_name = original_file_name
        self.plotting_data_evidence = trace_data
        self.ladder_plotting_data = ladder_data
        self.mapping_plotting_data = mapping_data
        self.truth_data = False
        self.accuracy = 0

    def set_truth_and_accuracy(self, truth_alleles, true_species):
        self.truth_data = truth_alleles
        self.true_species = true_species
        self.accuracy = accuracy.allele_accuracy(self.alleles_called,self.truth_data)