import FSAreader
import analysis
import os

import InputFileReaders


def processFSAfile(FSAfile, panel_info):

    print("****** Processing " + FSAfile + " **********")

    dye_to_channel_mapping, all_collected_data = FSAreader.readFSAFile(FSAfile)

    output_dir="./tmp/"
    run_name =all_collected_data["SpNm1"]
    run_folder = os.path.join(output_dir,run_name)
    if not(os.path.exists(run_folder )):
            os.makedirs(run_folder )


    relevant_loci = InputFileReaders.figure_out_loci_from_run_name(panel_info, run_name)

    if (relevant_loci == "FAIL" ):
        print("Uh-oh!  Can't figure out what panel to use for this FSA file!!")
    else:
        print("relevant_loci=", str(relevant_loci ))


    sixteen_peaks, threshold = analysis.getLadderPeaks(run_folder, all_collected_data)

    mappingFxn, left_domain_limit, right_domain_limit = \
        analysis.buildInterpolationdBasedOnLadder(run_folder, sixteen_peaks)

    analysis.RemapLadder(run_folder, all_collected_data ,mappingFxn,
                     left_domain_limit, right_domain_limit,
                     sixteen_peaks, threshold )


    # Channels we care about are ones with dyes in our panel.
    channels = set([loci_info_dict["dye"] for loci_info_dict in relevant_loci.values()])

    for channel in channels:
        analysis.RemapDataTrace(run_folder,
                        relevant_loci, #ie, the loci for this primer set
                        all_collected_data ,mappingFxn,
                        left_domain_limit, right_domain_limit,
                        sixteen_peaks, dye_to_channel_mapping[channel] )

