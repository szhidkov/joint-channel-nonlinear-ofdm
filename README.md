# Joint Channel Estimation and Nonlinear Distortion Compensation in OFDM Receivers 

This is a communication system simulation code for joint channel estimation and nonlinear distortion compensation in OFDM receivers (Bussgang-type receiver). For more details, see Ref [1].

## Getting Started

This simulation code is written in MATLAB/Octave.

### Prerequisites
NB: you might need vitdec() function from Communications Toolbox to run this simulation example.

## Running the simulation

The main simulation loop script is "main_sim.m". You need to specify your simulation parameters directly in the script. 

### Simulation parameters

See comments in "main_sim.m". 

### Getting results

The simulation results are plotted interactively during simulation and saved into "results.mat" after simulation has finished.

## References

[1] S. V. Zhidkov, "Joint Channel Estimation and Nonlinear Distortion Compensation in OFDM Receivers," Arxiv, 2016
https://arxiv.org/pdf/1612.09222v1  (see also presentation at: [https://cifrasoft.com/people/szhidkov/papers/zhidkov_rx_nonlinearity_compensation_ofdm.pdf](https://cifrasoft.com/people/szhidkov/papers/zhidkov_rx_nonlinearity_compensation_ofdm.pdf)


## Authors

* Sergey Zhidkov - *Initial work* ([https://cifrasoft.com/people/szhidkov/](https://cifrasoft.com/people/szhidkov/))


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


