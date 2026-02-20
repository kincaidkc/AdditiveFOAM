import re
import matplotlib.pyplot as plt
import numpy as np

#with plt.style.context(['science', 'nature']):

def parse_log_file(log_file_path):
    """
    Parses a log file to extract iteration numbers and average load imbalance values.

    Args:
        log_file_path (str): The path to the log file.

    Returns:
        tuple: A tuple containing two lists: iterations and load_imbalances.
    """
    iterations = []
    load_imbalances = []
    iteration_counter = 1

    log_pattern = re.compile(r"Average load imbalance: (\d+\.\d+)")
    
    try:
        with open(log_file_path, 'r') as f:
            for line in f:
                match = log_pattern.search(line)
                if match:
                    load_imbalances.append(float(match.group(1)))
                    iterations.append(iteration_counter)
                    iteration_counter += 1
    except FileNotFoundError:
        print(f"Error: The file '{log_file_path}' was not found.")
        return [], []

    return iterations, load_imbalances

def plot_load_imbalance(iterations, load_imbalances):
    """
    Plots the load imbalance as a function of iteration.

    Args:
        iterations (list): A list of iteration numbers.
        load_imbalances (list): A list of load imbalance values.
    """
    if not iterations or not load_imbalances:
        print("No data to plot.")
        return

    fig, ax = plt.subplots()
    ax.plot(iterations, load_imbalances, marker='o', linestyle='-')
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Average Load Imbalance")
    ax.set_title("Load Imbalance Over Iterations")
    ax.autoscale(tight=True)
    fig.savefig('load_imbalance_plot.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    log_file = 'log.additiveFoam'

    iterations, load_imbalances = parse_log_file(log_file)

    plot_load_imbalance(iterations, load_imbalances)
