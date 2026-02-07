# QoS-Aware Anti-Jamming Communication in Fluid Antenna Systems: Continuous or Discrete Position Design?

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Status: Research-in-Progress](https://img.shields.io/badge/Status-Research--in--Progress-orange.svg)]()

This repository contains the core algorithm implementations for the paper **"QoS-Aware Anti-Jamming Communication in Fluid Antenna Systems: Continuous or Discrete Position Design?"**. 

This research investigates the potential of **Fluid Antenna Systems (FAS)** to enhance anti-jamming performance in Multi-User MIMO (MU-MIMO) networks by leveraging the new degrees of freedom (DoF) provided by positional flexibility.

[Image of fluid antenna system architecture for anti-jamming communication]

---

## 🌟 Research Highlights

- **Framework**: A novel QoS-aware anti-jamming communication framework for MU-MIMO networks.
- **Robustness**: Formulates a sum-rate maximization problem considering **imperfect Jammer CSI** and practical antenna constraints.
- **Dual Design Perspectives**:
    - **Continuous Position Design**: Solved via an alternating optimization (AO) framework integrating Dinkelbach’s method, SCA, and Majorization-Maximization (MM).
    - **Discrete Position Design**: Reformulated into a regularized least squares form based on MMSE and solved using Block Coordinate Descent (BCD) and **Simultaneous Orthogonal Matching Pursuit (SOMP)**.
- **Performance**: Numerical results demonstrate that FAS significantly outperforms conventional Fixed-Position Antenna (FPA) systems in hostile jamming environments.

---

## 📂 Repository Structure

| Module | Description | Key Algorithms |
| :--- | :--- | :--- |
| `FAS_Continuous/` | Optimization for continuous antenna movement. | AO, Dinkelbach, SCA, MM |
| `FAS_Discrete/` | Optimization for selection from discrete ports. | MMSE, BCD, SOMP |
| `Core_Utils/` | Common utilities for channel modeling and jamming scenarios. | Robust CSI Processing |

---

## 🚀 Getting Started

### Prerequisites
- Python 3.8+
- NumPy, SciPy
- CVXPY (optional, for specific optimization benchmarks)

### Usage
Currently, the core logic for the proposed algorithms is provided in the `src/` directory. 
```bash
# Example: Run the continuous position optimization
python src/main_continuous_opt.py
