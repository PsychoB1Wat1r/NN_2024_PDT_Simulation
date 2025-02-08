## NFOFC for TDNNS with PDT switching: A system mode and time scheduler dual-dependent design

This repository contains the code implementation of our NN_2024 paper [here](https://www.sciencedirect.com/science/article/abs/pii/S0893608023006305).

If you find this repository useful, please cite our paper.
```
  @article{ma2023asynchronous,
  title={Non-fragile output-feedback control for time-delay neural networks with persistent dwell time switching: A system mode and time scheduler dual-dependent design},
  author={Zhou, Jianping and Ma, Xiaofeng and Yan, Zhilian and Arik, Sabri},
  journal={Neural Networks},
  volume={169},
  pages={733-743},
  year={2024},
  publisher={Elsevier}}
```
## Explanation:
We use **MATLAB** as the programming language and employ the **Mosek solver** along with the **YALMIP toolbox** as development tools to perform numerical solution and simulation for the proposed **system mode and time scheduler dual-dependent design** method.

## Mosek Solver:
Mosek is a high-performance mathematical optimization solver specifically designed for convex optimization problems. For official documentation, please refer to this [link](https://www.mosek.com/documentation/).

## YALMIP Toolbox:
YALMIP (Yet Another LMI Parser) is an optimization modeling toolbox for MATLAB that facilitates the formulation and solution of various optimization problems. Developed by [**Johan Löfberg**](https://scholar.google.com/citations?user=No-9sDUAAAAJ&hl=en), YALMIP is primarily used for convex optimization but can also handle certain non-convex problems.

It is important to note that YALMIP itself is **not a solver**; rather, it serves as a modeling interface that translates user-defined optimization problems into a standard form and then calls external solvers such as **Mosek**, **Gurobi**, **SDPT3**, **SeDuMi**, and others for solution. For official documentation, please refer to this [link](https://yalmip.github.io/).
