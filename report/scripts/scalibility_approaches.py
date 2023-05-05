import numpy as np
import matplotlib.pyplot as plt

# Data
# sphere3_res10_sequential_runtimes = np.array([1.1, 1.1, 1.1, 1.1, 1.1])
sphere10_res10_sequential_runtimes = np.array([44212, 47409, 42995])
sphere100_res10_sequential_runtimes = np.array([4.26949e+6, 4.05528e+6])
sphere200_res10_sequential_runtimes = np.array([2.31712e+8, 2.12297e+8])
sphere1000_res10_sequential_runtimes = np.array([0])
# sphere3_res10_parallel_runtimes = np.array([595, 1150, 1456, 1103, 1072, 982, 1077, 1138])
sphere10_res10_parallel_runtimes = np.array([1141, 1097, 1339, 981])
sphere100_res10_parallel_runtimes = np.array([12224, 12509, 11985, 12161])
sphere200_res10_parallel_runtimes = np.array([623369, 547624])
sphere1000_res10_parallel_runtimes = np.array([1.05086e+6, 1.06331e+6, 1.07816e+6, 1.09698e+6, 1.1033e+6])

# Process Data
# sphere3_res10_runtime = np.mean(sphere3_res10_sequential_runtimes)
sphere10_res10_sequential_runtime = np.mean(sphere10_res10_sequential_runtimes)
sphere100_res10_sequential_runtime = np.mean(sphere100_res10_sequential_runtimes)
sphere200_res10_sequential_runtime = np.mean(sphere200_res10_sequential_runtimes)
sphere1000_res10_sequential_runtime = np.mean(sphere1000_res10_sequential_runtimes)
# sphere3_res10_parallel_runtime = np.mean(sphere3_res10_parallel_runtimes)
sphere10_res10_parallel_runtime = np.mean(sphere10_res10_parallel_runtimes)
sphere100_res10_parallel_runtime = np.mean(sphere100_res10_parallel_runtimes)
sphere200_res10_parallel_runtime = np.mean(sphere200_res10_parallel_runtimes)
sphere1000_res10_parallel_runtime = np.mean(sphere1000_res10_parallel_runtimes)
sphere10_speedup = sphere10_res10_sequential_runtime / sphere10_res10_parallel_runtime
sphere100_speedup = sphere100_res10_sequential_runtime / sphere100_res10_parallel_runtime
sphere200_speedup = sphere200_res10_sequential_runtime / sphere200_res10_parallel_runtime
sphere1000_speedup = sphere1000_res10_sequential_runtime / sphere1000_res10_parallel_runtime

# Plot
plt.figure(figsize=(8, 4))
X = [10, 100, 200, 1000]
Y_sequential = [sphere10_res10_sequential_runtime, sphere100_res10_sequential_runtime, sphere200_res10_sequential_runtime, sphere1000_res10_sequential_runtime]
Y_parallel = [sphere10_res10_parallel_runtime, sphere100_res10_parallel_runtime, sphere200_res10_parallel_runtime, sphere1000_res10_parallel_runtime]
Y_speedup = [sphere10_speedup, sphere100_speedup, sphere200_speedup]
X_axis = np.arange(len(X))
plt.yscale('log')
plt.ylabel('Runtime (s)')
plt.xlabel('Number of Meshes')
plt.xticks(X_axis, X)
plt.bar(X_axis - 0.2, Y_sequential, 0.4, label = 'Sequential')
plt.bar(X_axis + 0.2, Y_parallel, 0.4, label = 'Parallel')
plt.legend(loc='upper left')
plt.twinx()
plt.plot([0,1,2], Y_speedup, 'r', label = 'Speedup')
plt.ylabel('Speedup')
plt.legend(loc='upper right')
plt.title('Runtime and Speedup of Parallelized MCD \n for different number of meshes')
plt.savefig('scalibility_approaches.png')
plt.show()