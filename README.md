# Electromagnetic Modelling of Wireless Communication Channels Assisted with Reconfigurable Intelligent Surfaces (RIS)

This repository contains MATLAB code for the senior design project focused on electromagnetic modelling of wireless communication channels enhanced by Reconfigurable Intelligent Surfaces (RIS). The project was completed in June 2022 at the Istanbul Technical University, Electrical-Electronics Faculty, Electronics and Communication Engineering Department.

**Abstract:**

This project investigates the impact of RIS on wireless communication channels by developing a physical channel model valid for sub-6GHz and mmWave bands. The model incorporates both deterministic RIS responses and stochastic channel parameters, including path loss, shadowing, multipath fading, polarization, and antenna radiation patterns.  The electromagnetic response of the RIS is modeled using a transmission line model.  The model assumes a Single Input Single Output (SISO) communication link with stationary transmitter and receiver. Simulations analyze the achievable rate's dependency on RIS size, positioning, and transmit power, comparing the results with existing literature.  A constrained optimization problem is formulated to maximize achievable rate with given design parameters and constraints.  The results provide a realistic perspective on the potential and limitations of using RIS in wireless communication.

**Key Features:**

* **Physical Channel Model:** Combines RIS EM response with a 3D polarized mmWave clustered channel model.
* **Transmission Line Model:**  Used for modelling the electromagnetic response of individual RIS elements.
* **Polarization and Antenna Patterns:**  Incorporates polarization-dependent fading and antenna radiation patterns.
* **Stochastic Channel Parameters:**  Models realistic channel conditions with random parameters like path loss, shadowing, and multipath fading.
* **SISO Communication:**  Focuses on single-antenna transmitter and receiver.
* **Achievable Rate Analysis:**  Simulations demonstrate the impact of RIS on achievable rate.
* **Optimization Problem:**  Formulates a constrained optimization problem for maximizing achievable rate.

**MATLAB Code:**

The repository includes MATLAB scripts and functions for:

* Generating stochastic channel parameters.
* Modelling the RIS response using the transmission line model.
* Calculating the overall channel coefficient.
* Simulating achievable rate under various scenarios.
* Formulating and solving the optimization problem.

**Simulation Results:**

The simulations explore the impact of:

* SNR on achievable rate.
* RIS size on achievable rate.
* RIS positioning on achievable rate.

**Future Work:**

* Solving the formulated optimization problem for RIS response.
* Validating the model with real-world measurements.
* Incorporating more sophisticated surface models (e.g., SEM-based models).
* Analyzing RIS response for different modulation schemes.
* Considering the mutual coupling between RIS elements.

**Further Information:**

For a more detailed description of the model and simulation results, please refer to the complete thesis document (available upon request).