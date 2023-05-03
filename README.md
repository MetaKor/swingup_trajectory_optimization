# swingup_trajectory_optimization
Optimal trajectory generation for the cartpole swing-up problem, formulated as a nonlinear program (NLP) in MATLAB.

Written by Luke Raus for Advanced Algorithms at Olin College.

## Details
To dig deeper, read the [full writeup PDF](https://github.com/MetaKor/swingup_trajectory_optimization/blob/main/Trajectory%20Generation%20via%20Optimization.pdf) or open the code! There are lots of comments :)

## Usage
This project runs in MATLAB and uses the `fmincon` NLP solver included in the Optimization Toolbox.

`example_swingup_optimization.m` contains a ready-to-run example problem for testing the implementation.

This script calls the `generate_swingup_trajectory` function to perform the actual optimization, and then validates this solution using the `simulate_cartpole` function.

## Results
The optimization returns open-loop control sequences and associated state trajectories as functions of time:

![Optimization results](/results/B_trajectory_N_50_m_03.png)

Enjoy!
