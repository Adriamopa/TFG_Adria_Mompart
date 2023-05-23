This whole pipeline has the usage of computing and analyzing iHS and std_iHS, outputting CSVs and plots, and returning the K-value and number of high values. All of this to consider different simulating or real conditions to compare iHS values and analyze the output. It can also change file format from VCF to ZARR.

The K-Value and high values are defined:
1. Values of standardized iHS greater than 2 (|iHS| >= 2) are considered as high values.
2. The number of high value SNPs are counted by sliding windows of 100kb of width and divided by the total number of SNPs in a given window.
3. The top 1% of windows are retrieved.
4. A new window is formed by the union of the windows that form the 1%.
5. The calculation of the number of high value SNPs divided by the number of SNPs is performed for the window, giving the value K.

An example of the execution of the pipeline is the following:

python3 PATH\TFG\Execution_interface.py -s std_ihs -vcf 'PATH\TFG\zarrs_tfg\NAME.zarr' -m 0.05 -w 10000 -mp 499999 -b 20

In this case, the calculation and analysis is only of the std_iHS of the NAME.zarr file, considering only variants with MAF over 0.05, 
windows for grouping results and plotting of 10000bp, the mutation is on position 499999, and standardize using 20 bins.

More than one file can be processed at the same time, like so:

python3 PATH\TFG_definitiu\Execution_interface.py -s std_ihs -vcf 'PATH\TFG_definitiu\zarrs_tfg\NAME.zarr','PATH\TFG_definitiu\zarrs_tfg\NAME2.zarr','PATH\TFG_definitiu\zarrs_tfg\NAME3.zarr' -m 0.05 -w 10000 -mp 499999 -b 20

As for now, all output files and plots produced are saved with predetermined names in predetermined folders, but this is to be optional in the future. This project is ambitioned to also be able to analyze cross_population extended haplotype homozygosity (XP_EHH), but
up to date only ihs is supported. Execution can be of 1 or of multiple files. Besides, it considers a chromosome of 1Mb fixed length. It should be however updated in the future.

Lastly, to change a file from vcf to zarr:

python3 PATH\TFG_definitiu\Execution_interface.py -vcf 'PATH\TFG_definitiu\vcfs_tfg\NAME.vcf' -z 'PATH\TFG_definitiu\zarrs_tfg\NAME.zarr'

Required programs to execute this pipeline include:
- SLiM simulating engine
- Python:
  - Scikit-allel
  - Matplotlib
  - Numpy
  - Time
  - Subprocess
