import matplotlib
import matplotlib.pyplot as plt
from io import BytesIO

# 准备数据
threads = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
speedups = [1.98, 1.65, 2.45, 2.50, 3.28, 3.42, 4.04, 4.35, 4.98, 5.26, 5.83, 6.11, 6.68, 7.04, 7.55]

# 绘制折线图
plt.figure(figsize=(6,4))  # 设置画布大小
plt.plot(threads, speedups, marker='o', linestyle='-', color='blue', label='Speedup')
plt.title('Speedup vs. Thread Count', fontsize=12)
plt.xlabel('Number of Threads', fontsize=10)
plt.ylabel('Speedup', fontsize=10)
plt.xticks(threads)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()

# 如果想保存成文件，请取消下面注释
plt.savefig('speedup_vs_threads.png', dpi=300)

plt.show()
