# QoS-Aware Anti-Jamming Communication in Fluid Antenna Systems: Continuous or Discrete Position Design?

[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: Research-in-Progress](https://img.shields.io/badge/Status-Core_Algorithms_Released-orange.svg)](#)

---

## 📖 简介 (Introduction)

由于无线信道的固有广播特性，无线通信日益容易受到恶意干扰攻击。本研究针对 **多用户多输入多输出 (MU-MIMO)** 网络，提出了一种新型的 **流体天线系统 (FAS)** 辅助的 **服务质量 (QoS) 感知抗干扰通信框架**。

利用流体天线的位置灵活性带来的新自由度，本框架旨在克服传统固定位置天线 (FPA) 的性能瓶颈，在复杂干扰环境下实现系统和速率的最优设计。

---

## 📡 系统模型 (System Model)


![System Model](https://youke.xn--y7xa690gmna.cn/s1/2026/02/07/69872ffff3b5f.webp)

---

## 🛠️ 核心算法 (Core Algorithms)

本项目针对两种不同的天线设计范式提供了完整的优化方案，并在 **MATLAB R2023b** 环境下完成验证：

| 设计方案 | 关键技术栈 | 核心逻辑 |
| :--- | :--- | :--- |
| **连续位置设计** | AO, SCA, MM, Dinkelbach | 通过交替优化框架精确定位天线在连续空间中的坐标 |
| **离散位置设计** | MMSE, BCD, SOMP | 基于最小均方误差准则，在离散端口集合中进行快速稀疏选择 |

---

## ⚠️ 重要说明 (Important Notice)

> [!IMPORTANT]
> **代码发布状态**：
> 1. **核心算法**：本仓库目前已给出论文中涉及的所有**核心算法函数代码**，运行**MainContinuous.m**和**MainDisceret.m**。
> 2. **完整仿真**：包含参数初始化、所有对比实验 (Benchmarks) 以及绘图脚本在内的**完整仿真代码**，将于论文正式发表后立即更新。

---

## 📊 数值仿真结论 (Numerical Results)

仿真结果表明，与传统的固定位置天线 (FPA) 以及现有的基准算法相比，所提框架具有以下优势：
* 显著提升了系统在干扰环境下的 **和速率 (Sum Rate)**。
* 在 **干扰者 CSI 不完美** 的实际场景中表现出极强的鲁棒性。
* 验证了流体天线在抗干扰通信中的巨大潜力。

---

## 📝 引用 (Citation)

如果您在研究中参考了本项目的算法或思路，请引用：

```bibtex
@article{YourName2026QoS,
  title={QoS-Aware Anti-Jamming Communication in Fluid Antenna Systems: Continuous or Discrete Position Design?},
  author={Yifan Guo and Co-authors},
  journal={Working Paper / Under Review},
  year={2026}
}
