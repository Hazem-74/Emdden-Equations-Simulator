# Emdden Equations Simulator

## Overview
This repository provides a Jupyter Notebook for numerically solving the Lane-Emden equation, a fundamental differential equation in astrophysics that models the structure of self-gravitating, spherically symmetric polytropic fluids, which serve as idealizations for stellar interiors. The project's goal is to calculate and visualize various key physical profiles—such as density, pressure, temperature, and luminosity—as a function of radius within these theoretical stellar models.

## Features
*   **Numerical Integration:** Implements a finite-difference method to solve the second-order Lane-Emden differential equation.
*   **Stellar Interior Profiles:** Generates comprehensive radial profiles for density ($\rho$), pressure ($P$), temperature ($T$), and luminosity ($L$) based on a user-defined polytropic index `n`.
*   **Physical Constant Derivation:** Automatically computes central density, pressure, temperature, the polytropic constant `K`, and the radial scaling factor `α` from the numerical solution.
*   **Data Visualization:** Utilizes `matplotlib` to plot the normalized potential ($\theta$) and the derived physical profiles against the scaled radius, offering clear insights into the stellar structure.
*   **Modularity:** Easily adaptable to different polytropic indices and numerical parameters for exploring a variety of stellar models.

## Files

### `Lane-Emdden.ipynb`
This Jupyter Notebook contains the core implementation for the numerical solution of the Lane-Emden equation and the subsequent calculation of stellar interior profiles.

**Purpose:**
The primary purpose of `Lane-Emdden.ipynb` is to construct a numerical model of a polytropic star's internal structure. It takes the Lane-Emden equation as its foundation, which describes the density distribution within a star under the assumption of a polytropic equation of state ($P = K\rho^{1+1/n}$). From this fundamental solution, the notebook proceeds to derive and visualize the radial dependencies of essential physical quantities: density, pressure, temperature, and luminosity.

**Methods and Logic:**
1.  **Equation Reformulation:** The notebook begins by transforming the Lane-Emden equation, `(1/ξ^2) d/dξ (ξ^2 dθ/dξ) = -θ^n`, into a form more suitable for numerical integration: `d^2θ/dξ^2 = -(2/ξ) dθ/dξ - θ^n`.
2.  **Numerical Integration (Finite Difference):** An iterative, forward-stepping finite-difference scheme is employed to solve this second-order ordinary differential equation.
    *   Starting from the stellar center at `ξ = 10^-6` (to avoid singularity at `ξ=0`), with initial conditions `θ(0)=1` and `dθ/dξ|_0=0`, the solution steps outwards in small increments (`Δξ`).
    *   At each step, the values of `dθ/dξ` and `θ` are updated using the discretized equations.
    *   The integration continues until `θ` becomes negative, indicating that the stellar surface (where `θ=0`) has been traversed.
3.  **Surface Parameter Determination:** To precisely locate the stellar surface (`ξ_1` where `θ=0`) and determine the derivative `dθ/dξ` at this point, linear interpolation is performed between the last two computed points (one with `θ > 0` and the subsequent one with `θ < 0`).
4.  **Physical Constant Calculation:** Utilizing the interpolated `ξ_1` and `dθ/dξ|_1`, the notebook calculates key physical constants characteristic of the stellar model, including:
    *   The dimensionless constants `a_n`, `c_n`, `b_n`.
    *   The central density ($\rho_c$), central pressure ($P_c$), and central temperature ($T_c$).
    *   The polytropic constant `K`.
    *   The radial scaling factor `α`, which relates the dimensionless `ξ` to the physical radius `r`.
5.  **Radial Profile Generation:** With all constants determined, the notebook then generates radial profiles for:
    *   Density: `ρ(r) = ρ_c θ^n`
    *   Pressure: `P(r) = P_c θ^(n+1)`
    *   Temperature: `T(r) = T_c θ^n`
    *   Luminosity: This is calculated by integrating the energy generation rate, specifically for the proton-proton (pp) chain (`ε_pp = 2.6 × 10^(-37) X^2 ρ T^4`), over the stellar volume.
6.  **Visualization:** The results are graphically presented using `matplotlib`, showing the behavior of `θ` versus `ξ`, and `ρ`, `P`, `T`, and `L` as a function of the scaled physical radius (`r/R_sun`).

## Requirements
To run this project, you will need:
*   Python 3.x
*   Jupyter Notebook
*   The following Python libraries:
    *   `numpy`
    *   `matplotlib`

## Installation
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/Hazem-74/Emdden-Equations-Simulator.git
    cd Emdden-Equations-Simulator
    ```
2.  **Create and activate a virtual environment (recommended):**
    ```bash
    python -m venv venv
    source venv/bin/activate # On Windows, use `venv\Scripts\activate`
    ```
3.  **Install the required Python packages:**
    ```bash
    pip install numpy matplotlib jupyter
    ```

## Usage
1.  **Launch Jupyter Notebook:**
    Navigate to the project directory in your terminal and run:
    ```bash
    jupyter notebook
    ```
2.  **Open the Notebook:**
    A web browser interface will open. From there, navigate to and click on `Lane-Emdden.ipynb` to open it.
3.  **Execute Cells:**
    Run the cells sequentially within the notebook. You can modify parameters such as the polytropic index `n` (default `3.3`) and numerical integration step `Δξ` in the initial code cells to explore different stellar models and refine accuracy.

## Validation and Comparison
The notebook is structured to facilitate validation by allowing users to set stellar parameters (e.g., solar mass `M_sun`, solar radius `R_sun`, mean molecular weight `μ`, hydrogen mass fraction `X`) that are relevant to standard stellar models. The generated radial profiles for density, pressure, temperature, and luminosity can then be compared qualitatively and, where applicable, quantitatively against published data or theoretical predictions for a Standard Solar Model (or other polytropic models with `n=3.3`).

## Contributing
Contributions to this project are welcome! If you have suggestions for enhancements, bug fixes, or new features, please feel free to open an issue or submit a pull request on the GitHub repository.

## License
This project is open-source and distributed under the [MIT License](LICENSE).