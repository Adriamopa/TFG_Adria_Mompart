import allel
import zarr
import argparse
from iHS_calculation import iHS

'''
This project is ambitioned to also be able to analyze cross_population extended haplotype homozygosity (XP_EHH), but
up to date only ihs is supported. Execution can be of 1 or of multiple files.
'''
FUNCTION_IHS = "ihs"
FUNCTION_STD_IHS = "std_ihs"
FUNCTION_XP_EHH = "xp_ehh"
FUNCTION_STD_XP_EHH = "std_xp_ehh"
FUNCTION_ALL = [FUNCTION_IHS, FUNCTION_STD_IHS, FUNCTION_XP_EHH, FUNCTION_STD_XP_EHH]

while True:
    if __name__ == '__main__':
        # total_start_time = time.perf_counter_ns()
        formatter = lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=10)
        parser = argparse.ArgumentParser(
            formatter_class=formatter,
            description="Computes and analyzes iHS and std_iHS, outputting CSVs and plots, and returning the K-value "
                        "and number of high values. It can also change file format from VCF to ZARR.\n"
        )
        parser.add_argument('-s', '--statistics',
                            help="Diversity and divergence statistics to calculate separated by comas.\n"
                                 "Possible values right now are:\n"
                                 "\t " + FUNCTION_IHS + ": For computing integrated haplotype score.\n"
                                 "\t " + FUNCTION_STD_IHS + ": For computing standardized iHS.\n", type=str)
        parser.add_argument('-z', '--zarr',
                            help="If argument -vcf is in VCF file format, it changes file format from VCF to ZARR and "
                                 "saves it in -z. \n",
                            type=str)
        parser.add_argument('-vcf', '--vcf_file_path',
                            help="ZARR file containing VCF data. Alternatively, VCF file to change to ZARR with -z "
                                 "argument\n", type=str)
        parser.add_argument('-m', '--min_maf', help="The minimum minor allele frequency to filter SNPs. Default is "
                                                    "0.05\n", type=float)
        parser.add_argument('-b', '--n_bins', help="Number of bins for standardization of iHS\n")
        parser.add_argument('-w', '--windows', help="Size of windows in base pairs to consider to group values for "
                                                    "the different plots. For K-value calculation and high value number"
                                                    ", a fixed 100000 is used.\n")
        parser.add_argument('-mp', '--mut_pos', help="Position of the selected mutation to plot. Defaults to the "
                                                     "central position of the chromosome\n")

        args = parser.parse_args()

        # Check the vcf_file argument: If it is None, exit the execution.
        vcf_file = None
        if args.vcf_file_path is None:
            raise ValueError("Error: Missing VCF file path. Check --help for more information.")
        # If zarr argument is not None, consider the vcf_file files to be in vcf file format and change them into zarr.
        # Once this is done, break the execution of the pipeline.
        if args.zarr is not None:
            vcf_files = args.vcf_file_path.split(',')
            zarrs = args.zarr.split(',')
            for i in range(len(vcf_files)):
                allel.vcf_to_zarr(vcf_files[i], zarrs[i], fields=('calldata/GT', 'variants/POS'))
            break

        # Check the statistics that we want to compute, if none is given, return an error.
        statistics = []
        if args.statistics is None:
            raise ValueError("Error: Provided statistics doesn't match with any of the available.\n"
                             "Check --help for more information.")
        # Check if the statistics are in the available ones.
        else:
            statistics = args.statistics.split(',')
            for stat in statistics:
                if stat not in FUNCTION_ALL:
                    raise ValueError("Error: Statistic " + stat + " doesn't match with any of the available.\n"
                                                                  "Check --help for more information.")

        # Save the min maf argument if given, else just go with the default.
        min_maf = 0.05
        if args.min_maf is not None:
            min_maf = args.min_maf

        # Save the bins argument if given, else just go with None (the default, and will convert into a reasonable
        # number depending on MAF and number of SNPs).
        n_bins = None
        if args.n_bins is not None:
            n_bins = int(args.n_bins)

        # Save the windows argument if given, else just go with None (the default, and will convert to 100000).
        windows = None
        if args.windows is not None:
            windows = int(args.windows)

        # Save the mutation position argument if given.
        mut_pos = None
        if args.mut_pos is not None:
            mut_pos = int(args.mut_pos)

        # Once all the argument have been saved, get all file paths and open them.
        vcf_files = args.vcf_file_path.split(',')
        for file in vcf_files:
            # Save the name of the file, which will be used for the name of the outputted files.
            name = file[file.index('iHS'):file.index('.zarr')]
            z_file = zarr.open(file, mode='r')
            # Extract all the positions and genotipes from the zarr file.
            positions = allel.SortedIndex(z_file['variants/POS'])
            genotypes = allel.GenotypeArray(z_file['calldata/GT'])
            # Print the number of positions/SNPs and the number of genotypes.
            print(genotypes.shape)
            # Get the haplotypes from the genotypes
            haplotypes = genotypes.to_haplotypes()
            # Get the allele count from the haplotypes, but only the biallelic ones
            ac = haplotypes.count_alleles(max_allele=1)
            # Determine the minimum MAF threshold considering the number of haplotypes
            min_min = int((int(haplotypes.shape[1]) * min_maf) - 1)
            print(min_min)
            # Get a filter of the biallelic variants that are segregating and have a minimum MAF over the threshold
            filtt = ac.is_segregating() & (ac.min(axis=1) > min_min)
            # Filter the haplotypes, positions and allele count
            h_seg = haplotypes[filtt]
            pos_seg = positions[filtt]
            ac_seg = ac.compress(filtt, axis=0)
            # With all the data obtained and extracted, execute the iHS/std_iHS
            iHS(h_seg, pos_seg, ac_seg, min_maf, statistics, name, windows, mut_pos, n_bins)
        break
