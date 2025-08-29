import matplotlib.pyplot as plt
import numpy as np
import re
import glob

def parse_benchmark_output(file_pattern):
    data = {}
    files = glob.glob(file_pattern)
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                match = re.match(r'Operation: (\w+), Algorithm: (\w+), Chunk Size: (\d+), Average Time: ([\d.]+) seconds, Standard Deviation: ([\d.]+) seconds', line)
                if match:
                    operation, algorithm, chunk_size, avg_time, stddev = match.groups()
                    chunk_size = int(chunk_size)
                    avg_time = float(avg_time)
                    stddev = float(stddev)
                    if operation not in data:
                        data[operation] = {}
                    if algorithm not in data[operation]:
                        data[operation][algorithm] = []
                    data[operation][algorithm].append((chunk_size, avg_time, stddev))
    return data

def plot_benchmark(data, operation):
    plt.figure(figsize=(10, 6))
    for algorithm, results in data[operation].items():
        results.sort()
        chunk_sizes, avg_times, stddevs = zip(*results)
        plt.errorbar(chunk_sizes, avg_times, yerr=stddevs, label=algorithm, capsize=5, marker='o')

    plt.xlabel('Chunk Size')
    plt.ylabel('Average Time (seconds)')
    plt.title(f'Benchmark Results for {operation.capitalize()} Operation')
    plt.legend()
    plt.grid(True)
    plt.savefig(f'{operation}_benchmark.png')
    plt.show()

def main():
    file_pattern = 'earlycoll_benchmark_*.out'
    data = parse_benchmark_output(file_pattern)
    
    for operation in data.keys():
        plot_benchmark(data, operation)

if __name__ == '__main__':
    main()
