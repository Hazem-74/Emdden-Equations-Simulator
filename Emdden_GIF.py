import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image
import io

# Constants
G = 6.67430e-11
M_sun = 1.9885e30
R_sun = 6.9634e8
mu = 0.6
mH = 1.6735575e-27
kB = 1.380649e-23
X = 0.7

def solve_lane_emden(n, delta_xi=1e-3, xi0=1e-6):
    xi = xi0
    theta = 1.0
    dtheta_dxi = 0.0
    xi_values = [xi]
    theta_values = [theta]
    dtheta_values = [dtheta_dxi]

    while True:
        dtheta_new = dtheta_dxi - delta_xi * (2 / xi * dtheta_dxi + theta ** n)
        theta_new = theta + delta_xi * dtheta_new
        xi_new = xi + delta_xi

        xi_values.append(xi_new)
        theta_values.append(theta_new)
        dtheta_values.append(dtheta_new)

        if theta_new < 0:
            break

        xi = xi_new
        theta = theta_new
        dtheta_dxi = dtheta_new

    xi_prev, xi_final = xi_values[-2], xi_values[-1]
    theta_prev, theta_final = theta_values[-2], theta_values[-1]
    dtheta_prev, dtheta_final = dtheta_values[-2], dtheta_values[-1]

    xi_1 = xi_prev - theta_prev * (xi_final - xi_prev) / (theta_final - theta_prev)
    dtheta_1 = dtheta_prev + (dtheta_final - dtheta_prev) * (xi_1 - xi_prev) / (xi_final - xi_prev)

    a_n = -xi_1 ** 3 / dtheta_1
    rho_c = a_n * 3 * M_sun / (4 * np.pi * R_sun ** 3)
    c_n = 1 / dtheta_1 ** 2
    P_c = c_n * G * M_sun ** 2 / R_sun ** 4
    b_n = 1 / ((n + 1) * xi_1 * (-dtheta_1))
    T_c = b_n * G * M_sun * mu * mH / (kB * R_sun)

    K = (4 * np.pi) ** (1 / n) * G / ((n + 1) * xi_1 ** ((n + 1) / n) * (-dtheta_1) ** ((1 - n) / n)) * \
        M_sun ** ((n - 1) / n) * R_sun ** ((3 - n) / n)
    alpha = np.sqrt((n + 1) * K * rho_c ** (1 / n - 1) / (4 * np.pi * G))

    r_values = np.array(xi_values[:-1]) * alpha
    theta_arr = np.array(theta_values[:-1])
    rho_profile = rho_c * theta_arr ** n
    P_profile = P_c * theta_arr ** (n + 1)
    T_profile = T_c * theta_arr ** n
    epsilon_pp = 2.6e-37 * X ** 2 * rho_profile * T_profile ** 4

    L_profile = []
    L_cumulative = 0.0
    for i in range(len(r_values) - 1):
        dr = r_values[i + 1] - r_values[i]
        dL = epsilon_pp[i] * 4 * np.pi * r_values[i] ** 2 * rho_profile[i] * dr
        L_cumulative += dL
        L_profile.append(L_cumulative)
    L_profile.append(L_cumulative)

    return r_values / R_sun, theta_arr, rho_profile, P_profile, T_profile, np.array(L_profile)

frames = []
n_values = np.round(np.arange(0.1, 5.0, 0.1), 2)

for n in n_values:
    try:
        r, theta, rho, P, T, L = solve_lane_emden(n=n)

        # Create a figure with custom grid specifications
        fig = plt.figure(figsize=(12, 10))
        gs = GridSpec(3, 2, height_ratios=[1, 1, 1], width_ratios=[1, 1], hspace=0.4, wspace=0.3)

        # Left column: 3 subplots with equal height (1/3 of the column)
        ax1 = fig.add_subplot(gs[0, 0])  # First row, left column
        ax2 = fig.add_subplot(gs[1, 0])  # Second row, left column
        ax3 = fig.add_subplot(gs[2, 0])  # Third row, left column

        # Right column: 2 subplots with equal height (1/2 of the column)
        ax4 = fig.add_subplot(gs[0:2, 1])  # First two rows, right column
        ax5 = fig.add_subplot(gs[2, 1])  # Third row, right column

        # Plot the data
        ax4.plot(r, theta, color='navy')
        ax4.set_title(rf"$\theta(\xi)$  at n = {n}")
        ax4.set_ylabel(r"$\theta$")
        ax4.set_xlabel(r"$\xi$")
        ax4.grid(True)

        ax2.plot(r, rho, color='darkgreen')
        ax2.set_title(rf"$\rho(r)$ at n = {n}")
        ax2.set_ylabel(r"$\rho(r)$")
        ax2.set_xlabel(r"$R_\odot$")
        ax2.grid(True)

        ax1.plot(r, P, color='firebrick')
        ax1.set_title(rf"$P(r)$ at n = {n}")
        ax1.set_ylabel(r"$P(r)$")
        ax1.set_xlabel(r"$R_\odot$")
        ax1.grid(True)

        ax3.plot(r, T, color='darkorange')
        ax3.set_title(rf"$T(r)$ at n = {n}")
        ax3.set_ylabel(r"$T(r)$")
        ax3.set_xlabel(r"$R_\odot$")
        ax3.grid(True)

        ax5.plot(r, L, color='purple')
        ax5.set_title(rf"$L(r)$ at n = {n}")
        ax5.set_ylabel(r"$L(r)$")
        ax5.set_xlabel(r"$R_\odot$")
        ax5.grid(True)

        # Add a main title
        fig.suptitle(f"Lane-Emden Profiles — Polytropic Index n = {n}", fontsize=16)

        # Save the figure to a buffer
        buf = io.BytesIO()
        plt.savefig(buf, format="png")
        buf.seek(0)
        frames.append(Image.open(buf))
        plt.close(fig)

    except Exception as e:
        print(f"Skipped n = {n} due to error: {e}")
# Save GIF
if frames:
    frames[0].save("lane_emden_profiles.gif",
                   save_all=True,
                   append_images=frames[1:],
                   duration=200,
                   loop=0)
    print("GIF saved as 'lane_emden_profiles.gif'")
else:
    print("No frames generated.")
