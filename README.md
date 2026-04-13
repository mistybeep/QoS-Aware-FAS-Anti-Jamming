# FAS-aided Robust Anti-Jamming Communications: Continuous and Discrete Positioning Designs

[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: Core Algorithms Released](https://img.shields.io/badge/Status-Core_Algorithms_Released-orange.svg)](#)

---

## 🌐 Language / 语言
> **Click below to switch between Chinese and English versions.** > **点击下方选项卡切换中英文版本。**

<details>
<summary><b>🇨🇳 中文版</b></summary>

### 📖 简介
由于无线信道的固有广播特性，无线通信日益容易受到恶意干扰攻击。本研究针对 **多用户多输入多输出 (MU-MIMO)** 网络，提出了一种新型的 **流体天线系统 (FAS)** 辅助的 **服务质量 (QoS) 感知抗干扰通信框架**。

利用流体天线的位置灵活性带来的新自由度，本框架旨在克服传统固定位置天线 (FPA) 的性能瓶颈，在复杂干扰环境下实现系统和速率的最优设计。

### 📡 系统模型
![System Model](https://youke.xn--y7xa690gmna.cn/s1/2026/02/07/69872ffff3b5f.webp)

### 🛠️ 核心算法
本项目针对两种不同的天线设计范式提供了完整的优化方案，并在 **MATLAB R2023b** 环境下完成验证：

| 设计方案 | 关键技术栈 | 核心逻辑 |
| :--- | :--- | :--- |
| **连续位置设计** | AO, SCA, MM, Dinkelbach | 通过交替优化框架精确定位天线在连续空间中的坐标 |
| **离散位置设计** | MMSE, BCD, SOMP | 基于最小均方误差准则，在离散端口集合中进行快速稀疏选择 |

### ⚠️ 重要说明
> [!IMPORTANT]
> **代码发布状态**：
> 1. **核心算法**：本仓库目前已给出论文中涉及的所有核心算法函数代码，请运行 `MainContinuous.m` 和 `MainDiscrete.m`。
> 2. **完整仿真**：包含参数初始化、所有对比实验 (Benchmarks) 以及绘图脚本在内的完整仿真代码，将于论文正式发表后立即更新。

### 📊 数值仿真结论
仿真结果表明，与传统的固定位置天线 (FPA) 以及现有的基准算法相比，所提框架具有以下优势：
* 显著提升了系统在干扰环境下的 **和速率 (Sum Rate)**。
* 在 **干扰者 CSI 不完美** 的实际场景中表现出极强的鲁棒性。
* 验证了流体天线在抗干扰通信中的巨大潜力。

</details>

<details open>
<summary><b>🇺🇸 English Version</b></summary>

### 📖 Introduction
Wireless communications are increasingly vulnerable to malicious jamming due to the inherent broadcast nature of the medium. This research proposes a novel **Quality of Service (QoS)-aware anti-jamming framework** for **Multi-User Multiple-Input Multiple-Output (MU-MIMO)** networks, leveraged by **Fluid Antenna Systems (FAS)**.

By exploiting the additional degrees of freedom (DoF) provided by the positional flexibility of fluid antennas, this framework aims to overcome the performance bottlenecks of traditional **Fixed-Position Antennas (FPA)** and optimize the system sum-rate in hostile interference environments.

### 📡 System Model
![System Model](https://youke.xn--y7xa690gmna.cn/s1/2026/02/07/69872ffff3b5f.webp)

### 🛠️ Core Algorithms
The project provides comprehensive optimization solutions for two distinct antenna design paradigms, verified in **MATLAB R2023b**:

| Design Paradigm | Key Tech Stack | Core Logic |
| :--- | :--- | :--- |
| **Continuous Position** | AO, SCA, MM, Dinkelbach | Precisely determines antenna coordinates in continuous space via an Alternating Optimization (AO) framework. |
| **Discrete Position** | MMSE, BCD, SOMP | Enables fast sparse selection among a finite set of ports based on the Minimum Mean Square Error (MMSE) criterion. |

### ⚠️ Important Notice
> [!IMPORTANT]
> **Code Release Status**:
> 1. **Core Algorithms**: This repository currently contains all core algorithmic functions mentioned in the paper. To test them, run `MainContinuous.m` and `MainDiscrete.m`.
> 2. **Full Simulation**: The complete codebase, including parameter initialization, all benchmark comparisons, and plotting scripts, will be released immediately upon the formal publication of the paper.

### 📊 Numerical Results
Simulations demonstrate that compared to traditional FPAs and existing benchmarks, the proposed framework offers:
* **Significant Sum-Rate Enhancement** in high-interference scenarios.
* **Strong Robustness** against imperfect Channel State Information (CSI) of the jammer.
* **Empirical Validation** of FAS potential in future secure communication layers.

</details>

---

## 📝 Citation
If you use these algorithms or ideas in your research, please cite:

```bibtex
@article{FASQOS,
  title={FAS-aided Robust Anti-Jamming Communications: Continuous and Discrete Positioning Designs},
  author={Yifan Guo},
  year={2026}
}
