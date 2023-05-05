import numpy as np
import matplotlib.pyplot as plt

# Data
naive_runtimes = np.array([316, 315, 319, 321, 389])
naive_omp_runtimes = np.array([321, 325, 328, 330, 345])
aabb_runtimes = np.array([2067, 2047, 2068, 2056, 2330])
aabb_omp_runtimes = np.array([2035, 1897, 1885, 1956, 1958])

# Process Data
naive_runtime = np.mean(naive_runtimes)
naive_omp_runtime = np.mean(naive_omp_runtimes)
aabb_runtime = np.mean(aabb_runtimes)
aabb_omp_runtime = np.mean(aabb_omp_runtimes)

# Plot
plt.figure(figsize=(8, 4))
plt.bar([1, 2, 3, 4], [naive_runtime, naive_omp_runtime, aabb_runtime, aabb_omp_runtime], color=['r', 'g', 'b', 'y'])
plt.text(1, naive_runtime + 1e1, str(int(naive_runtime)) + " us", ha='center', va='bottom')
plt.text(2, naive_omp_runtime + 1e1, str(int(naive_omp_runtime)) + " us", ha='center', va='bottom')
plt.text(3, aabb_runtime + 1e1, str(int(aabb_runtime)) + " us", ha='center', va='bottom')
plt.text(4, aabb_omp_runtime + 1e1, str(int(aabb_omp_runtime)) + " us", ha='center', va='bottom')
plt.xticks([1, 2, 3, 4], ['Naïve', 'Naïve OMP', 'AABB', 'AABB OMP'])
plt.xlabel('Approach')
plt.ylabel('Runtime (microseconds)')
plt.title('Runtime of Different Mesh Pair Generation Approaches \n on benchmark-random (1000 randomly generated spheres, collision margin = 0.01)')
plt.savefig('report/figs/runtime_pairgen.png')
# plt.show()