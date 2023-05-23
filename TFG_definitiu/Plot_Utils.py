import matplotlib.pyplot as plt
import numpy as np


def ihsPlot(score, sites, ylabel, path, absolute=False, mut_pos=None):
    '''
    Function to make and scatter plot of iHS values

    Parameter
    -------
    score : iHS values (unstandarized or standarized)
    sites: positions corresponding to the iHS values
    ylabel: Label for the y axis
    path: Output path in which to save the plot
    absolute: Variable to plot scores with its sign or absolute values.
    mut_pos: Position in which the mutation producing the signature is found.

    Return
    -------
    Scatter plot
    '''
    if mut_pos is None:
        mut_pos = 499999
    if absolute is True:
        score = np.abs(score)
    # Scale text size to be more readable
    plt.rcParams.update({'font.size': 15})
    plt.figure(figsize=(20, 5))
    # Make the plot
    plt.plot(sites, score, linestyle=' ', marker='.', mfc='none')
    # Add a vertical line in the mutation site
    plt.axvline(mut_pos,color='red',ymin=0.01,ymax=0.99)
    # Plot 2 horizontal lines that determine the threshold of high values
    plt.axhline(2,linestyle='dashed',color='orange')
    if absolute is False:
        plt.axhline(-2, linestyle='dashed', color='orange')
    plt.grid(axis='y')
    # Label the x and y axis
    plt.xlabel('Position (bp)')
    plt.ylabel(ylabel)
    # Save the figure and show the plot
    plt.savefig(path)
    #plt.show()


def barplot(y_label, positions, averages, width, path, mut_pos=None):
    '''
    A bar plot
    :param y_label: Label for the y axis
    :param positions: Positions referring to the middle point of the windows
    :param averages: Average iHS value inside each window
    :param width: Window width
    :param path: Path to save the plot
    :param mut_pos: Mutation position to plot the line
    :return: The bar plot
    '''
    if mut_pos is None:
        mut_pos = 499999
    plt.figure(figsize=(20, 5))
    plt.bar(x=positions, height=averages, width=width)
    plt.axvline(mut_pos, color='red', ymin=0.01, ymax=0.99)
    plt.grid(axis='y')
    plt.xlabel('Position (bp)')
    plt.ylabel(y_label)
    plt.savefig(path)
    #plt.show()


def plot_thr(y_label,centers, values, path, mut_pos=None):
    '''
    Plot a scatter plot and save the output. Also draw a vertival line for the mutation site.
    :param y_label: Label for the y axis
    :param centers: Middle position of the windows used to calculate the values
    :param values: Values to plot, proportion of high SNP values, usually
    :param path: Output path for the plot
    :param mut_pos: Mutation site to plot a vertical line
    :return: Scatter plot with centers as x values and values as y values
    '''
    if mut_pos is None:
        mut_pos = 499999
    plt.figure(figsize=(20, 5))
    plt.plot(centers,values,linestyle=' ', marker='.', mfc='none')
    plt.axvline(mut_pos,color='red',ymin=0.01,ymax=0.99)
    plt.grid(axis='y')
    plt.xlabel('Position (bp)')
    plt.ylabel(y_label)
    plt.savefig(path)
    #plt.show()


def plot_joint(y_label1, y_label2, plotx, ploty, barx, bary, width, path, mut_pos=None):
    '''
    Joint plot of a bar plot (args y_label1, barx, bary, width) and a scatter plot (y_label2, plotx, ploty) and a
    vertical line in the mutation site. Both plots use different y axis scaling but the same x axis scaling.
    :param y_label1: Label for the bar plot y axis
    :param y_label2: Label for the scatter plot y axis
    :param plotx: Position data for the scatter plot
    :param ploty: Value data for the scatter plot
    :param barx: Position data for the bar plot
    :param bary: Value data for the bar plot
    :param width: With value of the bars of the bar plot
    :param path: Path to save the output figure
    :param mut_pos: Mutation position to plot the vertical line
    :return: The joint plot of both bar and scatter plots.
    '''
    if mut_pos is None:
        mut_pos = 499999
    # Instantiate both plots and the figure that contains them
    fig = plt.figure(figsize=(20, 5))
    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", sharex=ax, frame_on=False)
    # Modify the bar plot
    ax.bar(x=barx, height=bary, width=width, color='lightblue')
    ax.axvline(mut_pos,color='red',ymin=0.01,ymax=0.99)
    ax.set_xlabel('Position (bp)')
    ax.tick_params(axis='y', colors='steelblue')
    ax.set_ylabel(y_label1,color='steelblue')
    ax.grid(axis='y', color='deepskyblue')
    # Modify the scatter plot
    ax2.plot(plotx,ploty,linestyle=' ', marker='.', mfc='none', color='blueviolet')
    ax2.yaxis.tick_right()
    ax2.set_ylabel(y_label2, color='blueviolet')
    ax2.yaxis.set_label_position('right')
    ax2.tick_params(axis='y', colors='blueviolet')
    ax2.grid(axis='y', color='slateblue')
    # Save and show the plot
    plt.savefig(path)
    #plt.show()