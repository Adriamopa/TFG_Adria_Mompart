import allel
import numpy as np
import time
import Data_Utils
import Plot_Utils


def iHS(h_seg, pos_seg, ac_seg, min_maf, statistics, name, windows, mut_pos, n_bins):
    '''
       Function to take several arguments and compute iHS value, save the score in csv and then create output plots and
       K-Value

       Return
       -------
       Nothing, it just saves plots and files produced by the execution of the function.
       '''
    # Print the name of the file currently working on.
    print('\n'+name+'\n')
    # Count iHS execution time
    IHS_start_time = time.perf_counter_ns()
    print('\nFor iHS:')
    # Execute iHS with scikit-allel package
    iHS = allel.ihs(h_seg, pos_seg, include_edges=False, min_ehh=0.3, min_maf=min_maf, use_threads=True)
    IHS_end_time = time.perf_counter_ns()
    # Print executing time
    print(f"IHS execution time [seconds]: {(IHS_end_time - IHS_start_time) * 1e-9}")
    # If the ihs statistic is given as an argument
    if 'ihs' in statistics:
        # Save the filtered (no nan values) iHS statistic in a csv file in plot folder
        flt = np.isnan(iHS)
        masked = np.column_stack((iHS[np.logical_not(flt)], pos_seg[np.logical_not(flt)]))
        file = 'Plot_Folder/Csvs/'+name+'.csv'
        np.savetxt(file, masked, delimiter=',')
        # Create a scatter plot of the iHS scores for every SNPs and its position, and save it in a determined folder
        Plot_Utils.ihsPlot(iHS[np.logical_not(flt)], pos_seg[np.logical_not(flt)], ylabel="Unstandardized iHS",
                           path="Plot_Folder/iHS_plots/"+name+".png", mut_pos=mut_pos)
        # Open the saved csv file to use the values for more plotting
        values, positions = Data_Utils.open_csv(file)
        # This function computes the average |iHS| value in set 'windows' width windows, and plots the results in a bar
        # plot. It also returns the values for the window middle position and the average value for the window.
        middles, averages = Data_Utils.finestres(pos=positions, vals=values, size=windows, ylabel="Unstandardized iHS",
                                                 abso=True, path='Plot_Folder/Average_iHS_windows_plots/' + name + '_' +
                                                                 str(windows) + '_ww_column_plot.png', mut_pos=mut_pos)
        # Open the csv file again and save the values again
        values, positions = Data_Utils.open_csv(file)
        # Create a threshold and window width variables to be used. Threshold is the consideration for high values and
        # the window width is the width for the calculation.
        threshold = '2'
        window_width = '100000'
        # Calculate the number or proportion of high |iHS| values (over 2) for the different 100000bp windows
        centers, proportions, starts, ends, non_zero, totals = Data_Utils.valors_alts(pos=positions, vals=values,
                                                 size=int(window_width), ylabel='|iHS|>=' + threshold + ', ' +
                                                                                window_width + 'bp windows',
                                                 path='Plot_Folder/Proportion_high_values_plots/' + name + '_threshold_'
                                                      + threshold + '_plot.png', threshold=threshold, mut_pos=mut_pos)
        # Plot together the average |iHS| in windows and the proportion of high |iHS| values
        Plot_Utils.plot_joint(y_label1='Average |iHS|', y_label2='|iHS|>='+threshold, plotx=centers, ploty=proportions,
                              barx=middles, bary=averages, width=windows, path='Plot_Folder/Joint_plots/'+name+
                                                                               '_joint_plot.png')
    if 'std_ihs' in statistics:
        # Count and print the standardization time
        std_IHS_start_time = time.perf_counter_ns()
        print('\nFor standardized iHS:')
        std_IHS, bins = allel.standardize_by_allele_count(score=iHS, aac=ac_seg[:, 1], n_bins=n_bins)
        std_IHS_end_time = time.perf_counter_ns()
        print(f"Standardized IHS execution time [seconds]: {(std_IHS_end_time - std_IHS_start_time) * 1e-9}")
        # Save the file in its corresponding folder
        file = 'Plot_Folder/Csvs/std_'+name+'_'+str(n_bins)+'_bins.csv'
        stacked_positions = np.column_stack((std_IHS, pos_seg))
        np.savetxt(file, stacked_positions, delimiter=',')
        # Create a scatter plot of the std_iHS for every SNPs and its position, and save it in a determined folder
        Plot_Utils.ihsPlot(std_IHS, pos_seg, ylabel="Standardized iHS", path="Plot_Folder/iHS_plots/std_"+name+'_' +
                                                                             str(n_bins)+"_bins.png", mut_pos=mut_pos)
        # Open the saved csv file to use the values for more plotting
        values, positions = Data_Utils.open_csv(file)
        # This function computes the average |std_iHS| value in set 'windows' width windows, and plots the results in a
        # bar plot. It also returns the values for the window middle position and the average value for the window.
        middles, averages = Data_Utils.finestres(pos=positions, vals=values, size=windows, ylabel="Standardized iHS",
                                                 abso=True, path='Plot_Folder/Average_iHS_windows_plots/std_'+name+'_' +
                                                 str(windows) + '_ww_' + str(n_bins)+'_bins_column_plot.png',
                                                 mut_pos=mut_pos)
        # Open the csv file again and save the values again
        values, positions = Data_Utils.open_csv(file)
        # Create a threshold and window width variables to be used. Threshold is the consideration for high values and
        # the window width is the width for the calculation.
        threshold = '2'
        vals_fins = '100000'
        # Calculate the number or proportion of high |std_iHS| values (over 2) for the different 100000bp windows
        centers, proportions, starts, ends, non_zero, totals = Data_Utils.valors_alts(pos=positions, vals=values,
                                                size=int(vals_fins), ylabel='|std_iHS|>=' + threshold + ', ' + vals_fins
                                                + 'bp windows', path='Plot_Folder/Proportion_high_values_plots/std_' +
                                                name + '_threshold_'+threshold+'_'+str(n_bins)+'_bins_plot.png',
                                                threshold=threshold, mut_pos=mut_pos)
        # Plot together the average |std_iHS| in windows and the proportion of high |iHS| values
        Plot_Utils.plot_joint(y_label1='Average |std_iHS|', y_label2='|std_iHS|>=' + threshold, plotx=centers[non_zero],
                              ploty=proportions[non_zero], barx=middles, bary=averages, width=windows,
                              path='Plot_Folder/Joint_plots/std_'+name+'_joint_plot.png')
        # Open the csv file again and save the values again
        values, positions = Data_Utils.open_csv(file)
        # Now, use the finestres function to plot the number of high |std_iHS| values for every non-overlapping window
        # of vals_fins size
        middles_number, averages_number = Data_Utils.finestres(pos=positions, vals=values, size=int(vals_fins),
                                          ylabel="#|iHS|>=2", abso=True,
                                          path='Plot_Folder/Num_high_values_plots/std_'+name+'_'+str(n_bins) +
                                          'num_high_values_plot.png', mut_pos=mut_pos, threshold=2)
        # Calculate the total number of high values in the chromosome
        total_high_values = sum(averages_number)
        # Open the csv file again and save the values again
        values, positions = Data_Utils.open_csv(file)
        # Compute the union window of the 1% of windows with the highes proportion of high |std_iHS| values
        final_window, state = Data_Utils.join_windows(proportions, starts, ends, mut_pos)
        # Do the same for the number of high values instead of the proportion
        final_window_num, estat2 = Data_Utils.join_windows(totals, starts, ends, mut_pos)
        # Calculate the K-Value
        K_value = Data_Utils.K_value_calc(final_window, values, positions, threshold=threshold)
        # Calculate the pseudo K-value (raw number instead of proportions
        pseudo_K_value = Data_Utils.K_value_calc(final_window_num, values, positions, threshold=threshold, prop=False)
        # Open a csv and save the result of the k_value to store it with the rest of results.
        with open('Results/K_Value/K_values.csv', 'r+') as nun:
            f = nun.read()
            if name in f:
                return
            else:
                nun.writelines(['std_'+name+' , ', str(K_value)+' , ', str(final_window), state+'\n'])
        # Open a csv and save the result of the pseudo k_value to store it with the rest of results.
        with open('Results/High_Value_Number/high_value_number.csv', 'r+') as non:
            i = non.read()
            if name in i:
                return
            else:
                non.writelines(['std_'+name+' , ', str(pseudo_K_value)+' , ', str(total_high_values)+' , ', str(final_window_num), estat2+'\n'])
        print(K_value)
        print(pseudo_K_value)
