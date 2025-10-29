import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

# Physical constants (SI units)
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

# Initial solve
r, theta, rho, P, T, L = solve_lane_emden(n=3.3)

# Plot with subplots
fig, axes = plt.subplots(3, 2, figsize=(12, 12), sharex=True)
plt.subplots_adjust(left=0.1, bottom=0.15, hspace=0.4)

theta_line, = axes[0, 0].plot(r, theta, color='blue')
axes[0, 0].set_ylabel(r'$\theta$')
axes[0, 0].set_title("Lane-Emden Solution")
axes[0, 0].set_xlabel(r'Radius ($R_\odot$)')

rho_line, = axes[0, 1].plot(r, rho, color='green')
axes[0, 1].set_ylabel(r'$\rho$ (kg/m³)')
axes[0, 1].set_xlabel(r'Radius ($R_\odot$)')

P_line, = axes[1, 0].plot(r, P, color='purple')
axes[1, 0].set_ylabel('P (Pa)')
axes[1, 0].set_xlabel(r'Radius ($R_\odot$)')

T_line, = axes[1, 1].plot(r, T, color='orange')
axes[1, 1].set_ylabel('T (K)')
axes[1, 1].set_xlabel(r'Radius ($R_\odot$)')

L_line, = axes[2, 0].plot(r[:-1], L[:-1], color='red')
axes[2, 0].set_ylabel('L (W)')
axes[2, 0].set_xlabel(r'Radius ($R_\odot$)')

axes[2, 1].axis('off')  # Empty subplot

# Slider setup
ax_slider = plt.axes([0.2, 0.05, 0.65, 0.03])
n_slider = Slider(ax_slider, 'Polytropic Index n', 1.0, 5.0, valinit=3.3, valstep=0.1)

def update(val):
    n_val = n_slider.val
    r, theta, rho, P, T, L = solve_lane_emden(n=n_val)
    theta_line.set_data(r, theta)
    rho_line.set_data(r, rho)
    P_line.set_data(r, P)
    T_line.set_data(r, T)
    L_line.set_data(r[:-1], L[:-1])
    for ax_row in axes:
        for ax in ax_row:
            ax.relim()
            ax.autoscale_view()
    fig.canvas.draw_idle()

n_slider.on_changed(update)

plt.show()
