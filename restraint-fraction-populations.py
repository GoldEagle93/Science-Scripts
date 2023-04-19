#plots the timeseries of the active fraction values
#based on the analyze_remd analysis scripts

from meld import vault
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import progressbar
import itertools

bold = ['#fcff5d', '#7dfc00', '#0ec434', '#228c68', '#8ad8e8', '#235b54', '#29bdab', '#3998f5', '#37294f', '#277da7', '#3750db', '#f22020', '#991919', '#ffcba5', '#e68f66', '#c56133', '#96341c', '#632819', '#ffc413', '#f47a22', '#2f2aa0', '#b732cc', '#772b9d', '#f07cab', '#d30b94', '#edeff3', '#c3a5b4', '#946aa2', '#5d4c86', '#201923']
#plots the timeseries of the active fractions for a set of replicas
def visualize_parameters(store, num_params):
    n_steps = store.max_safe_frame-1
    noes = ['dummy-contacts-right.dat', 'dummy-contacts-left.dat']
    fig, axs = plt.subplots(2, figsize=(10,20), sharex=True)
    counter = 0
    for param,noe in zip(range(num_params),noes):
        params = get_parameters(store, n_steps, param)
        max_val = (np.amax(params))
        params = params/max_val
        # plt.figure(figsize=(10,10))
        for index in [0,4,9,14,19,24,29]: #sample of replicas to plot    
        #for index in range(30):
            y,binEdges=np.histogram(params[index],bins=100)
            bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
            axs[counter].plot(bincenters, y, label='replica '+ str(index), marker='.', color = bold[index], markersize=30, linestyle='None')
            # axs[counter].title(noes[counter])
            # plt.hist(params[index],bins=100, label='replica '+ str(index), histtype='barstacked')
        counter += 1
        # plt.legend(loc=0)
        plt.legend(bbox_to_anchor=(0, 1.05, 1, 0), loc="lower left", mode="expand", ncol=7)
        plt.xlabel('fraction')
        plt.ylabel('density')
    axs[0].title.set_text(noes[0])
    axs[1].title.set_text(noes[1])
    plt.savefig('restraint-fraction-populations.png', dpi=300, bbox_inches='tight')

#gets the parameter values from the store
def get_parameters(store, n_steps, paramID):
    bar = get_progress_bar('Loading data', n_steps).start()
    params = np.zeros((store.n_replicas,n_steps), dtype=float)
    for index, frame in enumerate(range(n_steps)):
        bar.update(index)
        params[:, index] = store.load_discrete_parameters(frame)[:, paramID]
    return params

#gests progress bar
def get_progress_bar(label, n_steps):
    widgets = [
        "{}: ".format(label),
        progressbar.Percentage(),
        " ",
        progressbar.Bar(),
        " ",
        progressbar.ETA(),
    ]
    bar = progressbar.ProgressBar(maxval=n_steps, widgets=widgets)
    return bar

store = vault.DataStore.load_data_store()
store.initialize(mode='r')
visualize_parameters(store, 2)
