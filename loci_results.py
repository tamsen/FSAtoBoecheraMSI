import accuracy
from signal_processing.elastic_ladder_analysis import Mapping_Info


class loci_results:

    def __init__(self, raw_alleles_called, filtered_calls_for_loci, final_alleles_called,
                 original_file_name, trace_data, ladder_data, mapping_data, ladder_status):
        self.raw_alleles_called = raw_alleles_called
        self.filtered_alleles_called = filtered_calls_for_loci
        self.final_alleles_called = final_alleles_called
        self.original_file_name = original_file_name
        self.plotting_data_evidence = trace_data
        self.ladder_plotting_data = ladder_data
        self.mapping_plotting_data = mapping_data
        self.ladder_status = ladder_status
        self.truth_data = False
        self.true_species = False
        self.raw_accuracy = 0
        self.final_accuracy = 0
        self.BMW_closest_sample_name = False
        self.BMW_closest_sample_species = False
        self.BMW_closest_sample_similarity_score = False

    def set_truth_and_accuracy(self, truth_alleles, true_species):
        self.truth_data = truth_alleles
        self.true_species = true_species
        self.raw_accuracy = accuracy.allele_accuracy(self.raw_alleles_called, self.truth_data)
        self.final_accuracy = accuracy.allele_accuracy(self.final_alleles_called, self.truth_data)

    def get_high_confidence_allele_calls(self):

        if self.ladder_status == Mapping_Info.LadderState.Good:
            return self.final_alleles_called
        else:
            return []

    def get_BMW_determination_string(self):

        #Entry = [closest_sample_name, closest_species, similarity_score]

        if self.BMW_closest_sample_species:
            return str(self.BMW_closest_sample_species) + \
                " (score: " + str(self.BMW_closest_sample_similarity_score) + ")"
        else:
            return False


    def set_BMW_determination(self, BMW_determination):

        #Entry = [closest_sample_name, closest_species, similarity_score]

        if BMW_determination:
            self.BMW_closest_sample_name = BMW_determination[0]
            self.BMW_closest_sample_species = BMW_determination[1]
            self.BMW_closest_sample_similarity_score = BMW_determination[2]

def get_final_allele_calls(loci_results: loci_results):
    return loci_results.final_alleles_called

def get_raw_allele_calls(loci_results: loci_results):
    return loci_results.raw_alleles_called

