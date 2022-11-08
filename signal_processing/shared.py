class peak_calling_parameters:

    def __init__(self, num_peaks_needed, kernel_size,
                 min_distance_between_peaks, min_peak_width, threshold_multiplier):
        self.num_peaks_needed = num_peaks_needed
        self.kernel_size = kernel_size
        self.min_distance_between_peaks = min_distance_between_peaks
        self.min_peak_width = min_peak_width
        self.threshold_multiplier = threshold_multiplier
