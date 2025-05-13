# To get video of sigma vs y
import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import multiprocessing as mp
from functools import partial
from pathlib import Path
import argparse
from matplotlib.ticker import StrMethodFormatter
import random

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

def gettingSigma(filename):
    exe = ["./getSigma", filename]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    lines = stderr.decode("utf-8").split("\n")
    ytemp, sigmatemp = [],[]
    
    for i in range(len(lines)):
        values = lines[i].split(" ")
        if values == ['']:
            pass
        else:
            ytemp.append(float(values[0]))
            sigmatemp.append(float(values[1]))

    Y = np.asarray(ytemp)
    SIGMA = np.asarray(sigmatemp)

    return Y, SIGMA


def process_timestep(ti, folder):    
    """Process a single timestep."""
    t = 0.01*ti
    snapshot_file = Path(f"intermediate/snapshot-{t:.4f}")
    output_file = folder / f"{int(t * 1000):08d}.png"

    if not snapshot_file.exists():
        print(f"{snapshot_file} not found!")
        return
    
    if output_file.exists():
        print(f"{output_file} already present!")
        return
    
    y, sigma = gettingSigma(snapshot_file)
    # Part to plot
    AxesLabel, TickLabel = [50, 35]
    fig, ax = plt.subplots()
    ax.set_title(f't = {t:.2f}')
    
    ax.set_ylim([0,1])
    ax.set_xlim([0,2])
    ax.plot(y, sigma, 'o', markeredgecolor='black')

    # plt.show()
    plt.savefig(output_file, bbox_inches="tight", dpi=250)
    plt.close()
    print(f"{ti+1} is done")

def main():
    parser = argparse.ArgumentParser(description="Process facets for bubbles in sheets.")
    args = parser.parse_args()
    
    nGFS = 10000

    folder = Path('Video_sigma')  # output folder

    if not folder.is_dir():
        os.makedirs(folder)
    
    # Prepare the partial function with fixed arguments
    process_func = partial(process_timestep, folder=folder)
    
    # Use all available CPU cores
    num_processes = 5 
    
    nGFS_list = list(range(nGFS))
    random.shuffle(nGFS_list)
    
    # Create a pool of worker processes
    with mp.Pool(processes=num_processes) as pool:
        # Map the process_func to all timesteps
        pool.map(process_func, nGFS_list)

if __name__ == "__main__":
    main()
