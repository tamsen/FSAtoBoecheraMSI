
class FSA_File_Results:
  def __init__(self, MSI_loci_results_by_loci):
    self.MSI_loci_results_by_loci = MSI_loci_results_by_loci


class MSI_loci_results:
    def __init__(self, alleles_called, plotting_data):
      self.alleles_called = alleles_called
      self.plotting_data_evidence = plotting_data