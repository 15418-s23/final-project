import numpy as np
import matplotlib.pyplot as plt

# Data
naive_runtimes = np.array([1.15333e+7, 1.08855e+7, 1.08131e+7, 1.09229e+7, 1.0769e+7])
naive_omp_runtimes = np.array([987348, 1.15971e+6, 1.07113e+6, 1.12603e+6, 1.06784e+6])
gjk_sequential_runtimes = np.array([15214, 14335, 9636, 9732, 13045])
gjk_sequential_hillclimb_runtimes = np.array([1849, 1641, 2286, 2337, 2022])
gjk_parallel_runtimes = np.array([941, 786, 704, 1279, 903])
gjk_parallel_batched_runtimes = np.array([623, 673, 777, 601, 672])

# Process Data
naive_runtime = np.mean(naive_runtimes)
naive_omp_runtime = np.mean(naive_omp_runtimes)
gjk_sequential_runtime = np.mean(gjk_sequential_runtimes)
gjk_sequential_hillclimb_runtime = np.mean(gjk_sequential_hillclimb_runtimes)
gjk_parallel_runtime = np.mean(gjk_parallel_runtimes)
gjk_parallel_batched_runtime = np.mean(gjk_parallel_batched_runtimes)

# Plot
plt.figure(figsize=(12, 8))
plt.bar([1, 2, 3, 4, 5, 6], [naive_runtime, naive_omp_runtime, gjk_sequential_runtime, gjk_sequential_hillclimb_runtime, gjk_parallel_runtime, gjk_parallel_batched_runtime], color=['r', 'r', 'g', 'g', 'b', 'b'])
plt.text(1, naive_runtime + 1e6, str(int(naive_runtime)) + " us", ha='center', va='bottom')
plt.text(2, naive_omp_runtime + 1e6, str(int(naive_omp_runtime)) + " us", ha='center', va='bottom')
plt.text(3, gjk_sequential_runtime + 1e3, str(int(gjk_sequential_runtime)) + " us", ha='center', va='bottom')
plt.text(4, gjk_sequential_hillclimb_runtime + 1e3, str(int(gjk_sequential_hillclimb_runtime)) + " us", ha='center', va='bottom')
plt.text(5, gjk_parallel_runtime + 1e2, str(int(gjk_parallel_runtime)) + " us", ha='center', va='bottom')
plt.text(6, gjk_parallel_batched_runtime + 1e2, str(int(gjk_parallel_batched_runtime)) + " us", ha='center', va='bottom')
plt.xticks([1, 2, 3, 4, 5, 6], ['Naïve', 'Naïve OMP', 'GJK Sequential', 'GJK Sequential Hillclimb', 'GJK Parallel', 'GJK Parallel Batched'])
plt.xlabel('Approach')
plt.ylabel('Runtime (microseconds)')
plt.yscale('log')


plt.twinx()
naive_speedup = naive_runtime / naive_runtime
naive_omp_speedup = naive_runtime / naive_omp_runtime
gjk_sequential_speedup = naive_runtime / gjk_sequential_runtime
gjk_sequential_hillclimb_speedup = naive_runtime / gjk_sequential_hillclimb_runtime
gjk_parallel_speedup = naive_runtime / gjk_parallel_runtime
gjk_parallel_batched_speedup = naive_runtime / gjk_parallel_batched_runtime
plt.plot([1, 2, 3, 4, 5, 6], [naive_speedup, naive_omp_speedup, gjk_sequential_speedup, gjk_sequential_hillclimb_speedup, gjk_parallel_speedup, gjk_parallel_batched_speedup], 'k--')
plt.text(1, naive_speedup + 1e-1, str(round(naive_speedup, 2)) + "x", ha='center', va='bottom')
plt.text(2, naive_omp_speedup + 1e-1, str(round(naive_omp_speedup, 2)) + "x", ha='center', va='bottom')
plt.text(3, gjk_sequential_speedup + 1e-1, str(round(gjk_sequential_speedup, 2)) + "x", ha='center', va='bottom')
plt.text(4, gjk_sequential_hillclimb_speedup + 1e-1, str(round(gjk_sequential_hillclimb_speedup, 2)) + "x", ha='center', va='bottom')
plt.text(5, gjk_parallel_speedup + 1e-1, str(round(gjk_parallel_speedup, 2)) + "x", ha='center', va='bottom')
plt.text(6, gjk_parallel_batched_speedup + 1e-1, str(round(gjk_parallel_batched_speedup, 2)) + "x", ha='center', va='bottom')
plt.ylabel('Speedup')
plt.yscale('log')

plt.title('Runtime and Speedup of Different MCD Approaches \n on benchmark-lite (two meshes, each with 2404 faces)')
plt.savefig('report/figs/runtime_approaches.png')
# plt.show()

