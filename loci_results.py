import accuracy


class loci_results:

    def __init__(self, raw_alleles_called, filtered_calls_for_loci, final_alleles_called,
                 original_file_name, trace_data, ladder_data, mapping_data):
        self.raw_alleles_called = raw_alleles_called
        self.filtered_alleles_called = filtered_calls_for_loci
        self.final_alleles_called = final_alleles_called
        self.original_file_name = original_file_name
        self.plotting_data_evidence = trace_data
        self.ladder_plotting_data = ladder_data
        self.mapping_plotting_data = mapping_data
        self.truth_data = False
        self.true_species = "to be determined"
        self.raw_accuracy = 0
        self.final_accuracy = 0
        self.BMW_determination = False

    def set_truth_and_accuracy(self, truth_alleles, true_species):
        self.truth_data = truth_alleles
        self.true_species = true_species
        self.raw_accuracy = accuracy.allele_accuracy(self.raw_alleles_called, self.truth_data)
        self.final_accuracy = accuracy.allele_accuracy(self.final_alleles_called, self.truth_data)

    def set_BMW_determination(self, BMW_determination_line):

        if BMW_determination_line:
            splat = BMW_determination_line.split("\t")
            if len(splat) > 5:
                splat = BMW_determination_line.split("\t")
                score = splat[0]
                species = splat[2]
                self.BMW_determination = species + " (score: " + str(score) + ")"
