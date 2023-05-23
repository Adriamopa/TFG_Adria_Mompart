import csv
import numpy as np
import Plot_Utils


def open_csv(arxiu):
    '''
    This function opens a csv file with just 2 columns as values, removes all nan values and returns the columns.
    :param arxiu: a 2-column csv
    :return: 2 1D arrays for the first and second column of the csv
    '''
    rows = []
    with open(arxiu, 'r') as arx:
        lines = csv.reader(arx)
        for row in lines:
            rows.append(row)
    values = []
    positions = []
    for element in rows:
        if element[0] != 'nan':
            values.append(float(element[0]))
            positions.append(int(float(element[1])))
    print('Number of SNPs:',len(positions))
    return values, positions


def finestres(pos, vals, ylabel, path, mut_pos, size=None, abso=False, threshold=0):
    '''
    This function computes the average iHS value in set 'windows' width windows, and plots the results in a bar
    plot. It also returns the values for the window middle position and the average value for the window. If a threshold
     is given, it computes the number of values higher than the threshold per window.
    :param pos: Positions array
    :param vals: Values (iHS) array corresponding to de positions
    :param ylabel: Y label for the bar plot
    :param path: Path to save the figure of the bar plot
    :param mut_pos: Selected mutation position to plot
    :param size: Window size to compute
    :param abso: Whether iHS values have to be considered by the sign or as absolute values
    :param threshold: Threshold to consider for high value SNPs.
    :return: The center value of the window, and the value for the window
    '''
    # Firstly the size of the windows is considered
    if size is None:
        size = 100000
    # Then the windows are constructed
    inici = pos[0]
    final = pos[-1]
    fins = []
    ultim = False
    for inici_fin in range(inici, final, int(size)):
        parada = inici_fin + size
        if parada >= final:
            # last window
            parada = final
            ultim = True
        else:
            parada -= 1

        fins.append([inici_fin, parada])

        if ultim:
            break
    values = []

    # Then for every window, the average iHS value for the SNPs inside it is computed and saved in mitges object
    for fin in fins:
        ini = fin[0]
        fi = fin[1]
        num_elements = 0
        suma_valors = 0
        while ini <= pos[0] <= fi:
            if abs(vals[0]) >= float(threshold):
                num_elements += 1
                if abso:
                    suma_valors += abs(vals[0])
                else:
                    suma_valors += vals[0]
            if pos[0] == final:
                break
            pos.pop(0)
            vals.pop(0)
        if (num_elements != 0) and (threshold > 0):
            values.append(num_elements)
        elif num_elements != 0:
            values.append(suma_valors / num_elements)
        else:
            values.append(0)
    # The middle value of every window is calculated and stored in promig
    middle = []
    for element in fins:
        middle.append((element[0]+element[1])/2)
    # A bar plot is constructed with the data
    Plot_Utils.barplot(y_label=ylabel, positions=middle, averages=values, width=size, path=path, mut_pos=mut_pos)
    # Returns the middle point of the windows and the values for said windows
    return middle, values


def valors_alts(pos, vals, ylabel, path, mut_pos=None, size=None, threshold=2.0):
    '''
    Determine either the number of high values or the
    :param pos: Positions to
    :param vals:
    :param ylabel:
    :param path:
    :param mut_pos:
    :param threshold:
    :param size:
    :return:
    '''
    # Size is assessed
    if size is None:
        size = 100000
    # Windows are created, with inicis as beggining points and finals als ending points, and at the same time the
    # calculation takes place. Then, the number of high value SNPs is saved in els object and the proportion of high
    # value SNPs is saved in tot_els object.
    last = pos[-1]
    els = []
    tot_els = []
    inicis = []
    finals = []
    a = len(pos)
    for vent in range(a):
        ini = pos[0]
        fin = pos[0]+size
        num_elements = 0
        num_elements_totals = 0
        hm = 0
        while pos[hm] <= fin:
            num_elements_totals += 1
            if abs(vals[hm]) >= float(threshold):
                num_elements += 1
            if pos[hm] == last:
                break
            hm += 1
        pos.pop(0)
        vals.pop(0)
        if num_elements != 0:
            els.append(num_elements)
            tot_els.append(num_elements/num_elements_totals)
            inicis.append(ini)
            finals.append(fin)
        else:
            els.append(0)
            tot_els.append(0)
            inicis.append(ini)
            finals.append(fin)
        if fin >= last:
            break
    # With begin and end points, the middle point of each window is calculated.
    centers = []
    for i in range(len(inicis)):
        centers.append((int(inicis[i])+int(finals[i]))/2)
    # Create a filter to filter all 0 values
    non_zero = (np.array(tot_els) != 0)
    # Plot the windows and values, only if they are not zero
    Plot_Utils.plot_thr(y_label=ylabel, centers=np.array(centers)[non_zero], values=np.array(tot_els)[non_zero],
                        path=path, mut_pos=mut_pos)
    # The return returns the middle points of the windows used for calculation, the proportion of high values, the
    # beginning points of the windows, the end points of the windows, the filter for non_zero values, and the number of
    # high values per window
    return np.array(centers), np.array(tot_els), np.array(inicis), np.array(finals), non_zero, els


def join_windows(values, starts, ends, mut_pos, threshold_percentage=1):
    '''
    With the values, and starts and ends of windows, returns the joint window of the 'threshold percentage' top windows
    in terms of values. Returns the state of the window, whether it comprises the mutation site and whether it comprises
     all the windows in threshold percentage or not.
    :param values: proportion of high values SNPs for the window
    :param starts: starts of the windows
    :param ends: end points of the windows
    :param mut_pos: Selected mutation site
    :param threshold_percentage: Percentage of top windows to consider
    :return: Returns the state of the window, and the start and end points of the window.
    '''
    # Get the integer number of windows that make the top perentage
    num_windows = len(starts)//100*threshold_percentage
    # Create an array with value, start and end points and sort it for values. Then consider only the number of windows
    # that comprise the top percentage.
    vals_fins = np.array(values)
    vals_fins = np.column_stack((vals_fins,np.array(starts)))
    vals_fins = np.column_stack((vals_fins,np.array(ends)))
    vals_sorted = vals_fins[vals_fins[:, 0].argsort()][::-1]
    vals_definitive = vals_sorted[:num_windows]
    print(vals_definitive)
    # Only with the start and end points, use the find_overlapping_intervals function to find the window that is the
    # union of the others.
    windows = vals_definitive[:, 1:]
    final_window, state = find_overlapping_intervals(windows, mut_pos)
    # If the window does not comprise the mutation site, give a '*' notice.
    if final_window[0] <= mut_pos <= final_window[1]:
        return final_window, ''+state
    else:
        return final_window, ' , *'+state


def K_value_calc(window, values, positions, threshold=2, prop=True):
    '''
    Obtains the K-value. It needs a set of iHS values and its positions, and then the window in which to calculate the
    K-value. If prop is false, it instead calculates only the number of high values in the given window
    :param window: Window in which to calculate the K-Value
    :param values: SNP iHS values
    :param positions: positions corresponding to the iHS values
    :param threshold: threshold for the consideration of high values
    :param prop: If true, calculate the K-value, if false calculate the number of high values
    :return: K-value or number of high values in window
    '''
    threshold = float(threshold)
    suma_mes_2 = 0
    comptador_vals = 0
    print(window)
    print(len(positions))
    for element in range(len(positions)):
        if window[1] >= positions[element] >= window[0]:
            comptador_vals += 1
            if abs(values[element]) >= threshold:
                suma_mes_2 += 1
    if not prop:
        return suma_mes_2
    k_value = suma_mes_2/comptador_vals
    return k_value


def find_overlapping_intervals(intervals, mut_pos):
    """
    Given a list of intervals, returns the union of overlapping intervals. Return the state of the return window,
    whether it comprises all the windows or not
    :param mut_pos: Mutation position to consider.
    :param intervals: 2d Array containing start and end point of the windows to find overlapping intervals.
    :return: A tuple representing the start and end points of the union of overlapping intervals, and the state of the
    window.
    """
    # Change value type of elements in the interval into int
    intervals = intervals.astype(int)

    # Sort the intervals in ascending order
    intervals = intervals[intervals[:, 0].argsort()]

    # Initialize the current interval with the first interval in the list
    current_interval = intervals[0]

    # Initialize the list of overlapping intervals with the current interval
    overlapping_intervals = [current_interval]
    interval_sum = [1]

    # Loop through the remaining intervals
    for interval in intervals[1:]:
        # If the current interval overlaps with the next interval
        if current_interval[1] >= interval[0]:
            # Update the end time of the current interval
            current_interval = (current_interval[0], max(current_interval[1], interval[1]))
            # Add the updated current interval to the list of overlapping intervals
            overlapping_intervals[-1] = current_interval
            interval_sum[-1] += 1
        else:
            # Otherwise, set the next interval as the new current interval
            current_interval = interval
            # Add the new current interval to the list of overlapping intervals
            overlapping_intervals.append(current_interval)
            interval_sum.append(1)
    # Return the interval that was built using the highest amount of windows, and if there is a draw, return that which
    # comprises the mutation site.
    final_window = overlapping_intervals[0]
    size = interval_sum[0]
    for element in range(len(overlapping_intervals)):
        if (interval_sum[element]) > size:
            final_window = overlapping_intervals[element]
        elif (interval_sum[element] == size) and (overlapping_intervals[element][0] <= mut_pos <=
                                                  overlapping_intervals[element][1]):
            final_window = overlapping_intervals[element]
    state = ''
    # If the overlapping interval was not built using all the windows given, return a ! message to consider.
    if size != len(intervals):
        state = ' , !'
    return final_window, state
