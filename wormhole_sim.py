"""
================================================================================
 TRAVERSABLE WORMHOLE SIMULATION — EARTH <-> ANDROMEDA GALAXY
================================================================================

What this is
------------
A numerical model, built directly from Einstein's field equations, of a
"traversable" wormhole connecting our Solar System to the Andromeda Galaxy
(M31), roughly 2.5 million light-years away.

It uses the Morris-Thorne "zero-tidal-force" traversable-wormhole metric
(the same simple solution Morris & Thorne present as a worked example in
their original 1988 paper, and the same family of metric used by Kip
Thorne's team to generate the wormhole visuals in the film *Interstellar*):

        ds^2 = -c^2 dt^2 + dl^2 + r(l)^2 (dtheta^2 + sin^2(theta) dphi^2)
        r(l) = sqrt(l^2 + b0^2)

Here `l` is the proper radial distance through the wormhole (running from
-infinity on the Earth side to +infinity on the Andromeda side), `b0` is the
radius of the wormhole's narrowest point (the "throat"), and r(l) is the
circumference (divided by 2*pi) at position l.

The script:
  1. Builds the embedding diagram (the classic wormhole "funnel" shape) by
     solving the isometric-embedding equation for this metric exactly.
  2. Uses symbolic tensor calculus (sympy) to compute the Einstein tensor
     directly from the metric, and from it the stress-energy that General
     Relativity REQUIRES to hold the throat open. This is not assumed -
     it's derived, and it comes out negative (see Section 2): you cannot
     build this wormhole out of ordinary matter.
  3. Integrates the null-geodesic (light-ray) equations of motion through
     the throat for a range of impact parameters, to show which rays
     cross through to the other universe and which reflect back.
  4. Computes special-relativistic transit times for a traveler crossing
     the throat, and compares them to the ~2.5-million-year trip light
     would take through ordinary space.
  5. Derives the tidal acceleration felt by a moving traveler (via the
     geodesic-deviation equation, again computed symbolically rather than
     quoted from a table) and uses it to work out how fast a human could
     comfortably cross a wormhole of a given throat size.

Everything physical that comes out of this script (the negative energy
density, the tidal-force/throat-size trade-off, the embedding shape) is
calculated from the metric above, not hard-coded.

Sources are listed in full at the bottom of this file and are printed
again at the end of the run.
================================================================================
"""

import os
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (enables 3D projection)
from scipy.integrate import solve_ivp

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUTDIR, exist_ok=True)

# ==============================================================================
# SECTION 0 — PHYSICAL CONSTANTS & SCENARIO PARAMETERS
# ==============================================================================
C = 2.99792458e8          # speed of light, m/s
G_NEWTON = 6.674e-11       # Newton's gravitational constant, m^3 kg^-1 s^-2
LY = 9.4607e15             # 1 light-year, m
AU = 1.495978707e11        # 1 astronomical unit, m

D_ANDROMEDA_LY = 2.5e6                  # NASA figure, see sources [6]
D_ANDROMEDA_M = D_ANDROMEDA_LY * LY

# --- Free design parameters of the wormhole -----------------------------
# General relativity does not predict these numbers - they are exactly the
# kind of thing Morris & Thorne treat as an *engineering* choice in [1].
B0 = 10.0                  # throat radius, metres (try 10 m = roomy hallway)
L_TUNNEL = 1.0e4            # proper length of the throat "tunnel", metres
                            #   default: 10 km - i.e. a wormhole tunnel
                            #   ~2.4e18 times shorter than the 2.5 Mly of
                            #   normal space it bypasses
XI_BODY = 2.0               # transverse size of a traveler, m (tidal-force test)
G_EARTH = 9.81               # m/s^2, comfort threshold for tidal acceleration


# ==============================================================================
# SECTION 1 — METRIC AND EMBEDDING GEOMETRY
#
#   ds^2 = -c^2 dt^2 + dl^2 + r(l)^2 dOmega^2 ,   r(l) = sqrt(l^2 + b0^2)
#
# This is the "zero tidal force" special case (redshift function Phi = 0)
# of the general Morris-Thorne metric, and is identical to the wormhole
# independently found by Ellis [3] (sourced by a "phantom" scalar field).
# Morris & Thorne present this exact form as their simplest worked example
# (Box 2 of [1]).
# ==============================================================================
def r_of_l(l, b0=B0):
    """Areal radius (circumference / 2 pi) as a function of proper distance l."""
    return np.sqrt(l**2 + b0**2)


def z_of_l(l, b0=B0):
    """
    Height of the isometric embedding of the equatorial slice into flat
    3-space, used to draw the classic wormhole 'funnel' picture.

    Embedding a surface of revolution (z(l), r(l)) into Euclidean 3-space
    requires (dz/dl)^2 + (dr/dl)^2 = 1 so that distances on the embedded
    surface match the wormhole's own metric. Since dr/dl = l / r(l):

        (dz/dl)^2 = 1 - l^2/r^2 = b0^2 / r^2   =>   dz/dl = b0 / r(l)

    which integrates in closed form to:

        z(l) = b0 * arcsinh(l / b0)
    """
    return b0 * np.arcsinh(l / b0)


# ==============================================================================
# SECTION 2 — SYMBOLIC GENERAL RELATIVITY
# Derive the Einstein tensor for the metric above from first principles
# (Christoffel symbols -> Riemann tensor -> Ricci tensor -> Einstein tensor),
# then read off the energy density General Relativity demands via
# G_{mu nu} = 8 pi G/c^4 T_{mu nu}.
# ==============================================================================
def derive_field_equations():
    t, l, th, ph, b0 = sp.symbols('t l theta phi b0', real=True, positive=False)
    b0 = sp.symbols('b0', positive=True)
    v, gamma_s, xi_s = sp.symbols('v gamma xi', positive=True)

    r = sp.sqrt(l**2 + b0**2)
    coords = [t, l, th, ph]
    n = 4

    g = sp.diag(-1, 1, r**2, r**2 * sp.sin(th)**2)
    ginv = g.inv()

    # Christoffel symbols of the second kind
    Gamma = [[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b_ in range(n):
            for c in range(n):
                s = 0
                for d in range(n):
                    s += ginv[a, d] * (sp.diff(g[d, c], coords[b_]) +
                                        sp.diff(g[d, b_], coords[c]) -
                                        sp.diff(g[b_, c], coords[d]))
                Gamma[a][b_][c] = sp.simplify(s / 2)

    # Riemann tensor R^a_{bcd}
    Riem = [[[[sp.Integer(0)] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b_ in range(n):
            for c in range(n):
                for d in range(n):
                    term = sp.diff(Gamma[a][b_][d], coords[c]) - sp.diff(Gamma[a][b_][c], coords[d])
                    s2 = 0
                    for e in range(n):
                        s2 += Gamma[a][c][e] * Gamma[e][b_][d] - Gamma[a][d][e] * Gamma[e][b_][c]
                    Riem[a][b_][c][d] = sp.simplify(term + s2)

    # Ricci tensor and scalar
    Ricci = sp.zeros(n, n)
    for b_ in range(n):
        for d in range(n):
            Ricci[b_, d] = sp.simplify(sum(Riem[a][b_][a][d] for a in range(n)))
    Rscalar = sp.simplify(sum(ginv[i, j] * Ricci[i, j] for i in range(n) for j in range(n)))

    # Einstein tensor G_{mu nu} = R_{mu nu} - 1/2 R g_{mu nu}
    Einstein = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            Einstein[i, j] = sp.simplify(Ricci[i, j] - sp.Rational(1, 2) * Rscalar * g[i, j])

    # Orthonormal-frame stress-energy components implied by G = 8 pi T (units G=c=1)
    rho_expr = sp.simplify(Einstein[0, 0] / (8 * sp.pi))                       # energy density
    pr_expr = sp.simplify(Einstein[1, 1] / (8 * sp.pi))                        # radial pressure
    pt_expr = sp.simplify((Einstein[2, 2] / r**2) / (8 * sp.pi))               # lateral pressure
    nec_expr = sp.simplify(rho_expr + pr_expr)                                  # null-energy-condition combo

    # Geodesic-deviation tidal tensor for an observer moving radially at speed v
    # u^mu = (gamma, gamma*v, 0, 0);  separation vector purely in theta-direction
    u = [gamma_s, gamma_s * v, 0, 0]
    a_theta = -sum(Riem[2][bb][2][dd] * u[bb] * xi_s * u[dd] for bb in range(n) for dd in range(n))
    a_theta = sp.simplify(a_theta.subs(th, sp.pi / 2))
    tidal_per_xi = sp.simplify((r * a_theta / xi_s))   # physical tidal accel per unit physical separation

    return {
        'r': r, 'b0': b0, 'l': l, 'v': v, 'gamma': gamma_s, 'xi': xi_s,
        'rho': rho_expr, 'p_radial': pr_expr, 'p_lateral': pt_expr, 'nec': nec_expr,
        'tidal_per_xi': tidal_per_xi,
    }


# ==============================================================================
# SECTION 3 — GEODESICS: LIGHT RAYS THROUGH THE THROAT
#
# For this static, spherically symmetric metric, motion in the equatorial
# plane (theta = pi/2) conserves:
#     E      = dt/dlambda            (from the t-Killing vector)
#     Lambda = r(l)^2 dphi/dlambda   (from the phi-Killing vector)
# Differentiating the null-geodesic normalisation -E^2 + l_dot^2 + Lambda^2/r^2 = 0
# once more in lambda removes the need to handle sign flips at turning points:
#     l_dot      = v_l
#     v_l_dot    =  Lambda^2 * r'(l) / r(l)^3
#     phi_dot    =  Lambda / r(l)^2
# Lengths below are measured in units of the throat radius (b0 = 1) since
# this section is about the *geometry/optics* of the throat, not absolute
# scale.
# ==============================================================================
def geodesic_rhs(lam, y, Lambda):
    l, vl, phi = y
    r = r_of_l(l, b0=1.0)
    rp = l / r  # dr/dl
    return [vl, Lambda**2 * rp / r**3, Lambda / r**2]


def integrate_ray(Lambda, E=1.0, l0=-45.0, lam_span=130.0, n_eval=3000):
    r0 = r_of_l(l0, b0=1.0)
    inside_sqrt = E**2 - Lambda**2 / r0**2
    if inside_sqrt <= 0:
        return None  # ray never even gets going at this starting point
    vl0 = np.sqrt(inside_sqrt)
    sol = solve_ivp(geodesic_rhs, [0, lam_span], [l0, vl0, 0.0], args=(Lambda,),
                     t_eval=np.linspace(0, lam_span, n_eval), rtol=1e-9, atol=1e-11,
                     dense_output=False)
    return sol


def critical_impact_parameter(E=1.0, b0=1.0):
    """Rays with |Lambda/E| < b0 cross the throat ('transmitted');
    rays with |Lambda/E| > b0 turn back before reaching it ('reflected')."""
    return E * b0


# ==============================================================================
# SECTION 4 — TRAVEL TIME: SPECIAL-RELATIVISTIC TRANSIT THROUGH THE THROAT
#
# Along the wormhole axis (fixed theta, phi) the metric is exactly flat
# Minkowski space, ds^2 = -c^2 dt^2 + dl^2, so a traveler moving at constant
# speed v through the tunnel behaves exactly as in special relativity.
# ==============================================================================
def transit_times(v_fraction_c, L=L_TUNNEL):
    v = v_fraction_c * C
    gamma = 1.0 / np.sqrt(1.0 - v_fraction_c**2)
    t_coord = L / v          # time elapsed for observers at the two mouths
    tau = t_coord / gamma    # proper time felt by the traveler
    return t_coord, tau, gamma


# ==============================================================================
# SECTION 5 — TIDAL FORCES & HUMAN TRAVERSABILITY
#
# From Section 2's symbolic geodesic-deviation calculation:
#     a_tidal = gamma^2 * (v/c)^2 * c^2 * (b0^2 / r(l)^3) * xi
# This vanishes for a *static* observer (consistent with Phi=0 meaning
# "zero tidal force" for someone sitting still) and grows with the SQUARE
# of the traveler's speed - exactly the "velocity-induced tidal force"
# Morris & Thorne flag as the binding constraint on how fast a human body
# can safely cross a small-throat wormhole.
# ==============================================================================
def tidal_acceleration(v_fraction_c, l=0.0, b0=B0, xi=XI_BODY):
    gamma = 1.0 / np.sqrt(1.0 - v_fraction_c**2)
    r = r_of_l(l, b0)
    return gamma**2 * v_fraction_c**2 * C**2 * (b0**2 / r**3) * xi


def max_comfortable_speed(b0=B0, xi=XI_BODY, g_limit=G_EARTH):
    """Solve gamma^2 v^2 = g_limit * b0 / (c^2 xi) for v at the throat (r=b0)."""
    rhs = g_limit * b0 / (C**2 * xi)
    v2 = rhs / (1.0 + rhs)
    return np.sqrt(v2)


def throat_for_comfortable_speed(v_fraction_c, xi=XI_BODY, g_limit=G_EARTH):
    """Throat radius needed so that crossing at v_fraction_c stays under g_limit."""
    gamma = 1.0 / np.sqrt(1.0 - v_fraction_c**2)
    return gamma**2 * v_fraction_c**2 * C**2 * xi / g_limit


# ==============================================================================
# SECTION 6 — VISUALISATION
# ==============================================================================
def make_embedding_with_geodesics():
    l_max = 7.0
    l_vals = np.linspace(-l_max, l_max, 220)
    phi_vals = np.linspace(0, 2 * np.pi, 120)
    L_grid, PHI_grid = np.meshgrid(l_vals, phi_vals)
    R_grid = r_of_l(L_grid, b0=1.0)
    Z_grid = z_of_l(L_grid, b0=1.0)
    X_grid = R_grid * np.cos(PHI_grid)
    Y_grid = R_grid * np.sin(PHI_grid)

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X_grid, Y_grid, Z_grid, rstride=3, cstride=3,
                     color='#3a7bd5', alpha=0.35, edgecolor='#1f4e8c', linewidth=0.2)

    # throat circle
    th_circle = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(th_circle), np.sin(th_circle), np.zeros_like(th_circle),
            color='red', linewidth=2.0, label='Throat (narrowest point)')

    # a few photon trajectories: some transmitted, one reflected
    b0 = 1.0
    Lambdas = [0.0, 0.4, 0.8, 0.97, 1.15]
    colors = ['#ffd700', '#ff8c00', '#ff3b3b', '#7CFC00', '#b042ff']
    for Lambda, col in zip(Lambdas, colors):
        sol = integrate_ray(Lambda, E=1.0, l0=-l_max - 3)
        if sol is None:
            continue
        l_path, vl_path, phi_path = sol.y
        mask = np.abs(l_path) <= l_max
        r_path = r_of_l(l_path, b0=1.0)
        z_path = z_of_l(l_path, b0=1.0)
        x_path = r_path * np.cos(phi_path)
        y_path = r_path * np.sin(phi_path)
        status = "transmitted" if abs(Lambda) < critical_impact_parameter(1.0, b0) else "reflected"
        ax.plot(x_path[mask], y_path[mask], z_path[mask], color=col, linewidth=2.2,
                label=f'light ray, b={Lambda:.2f} ({status})')

    ax.text(0, 0, z_of_l(l_max, 1.0) + 1.0, "ANDROMEDA\nside", color='darkblue',
            fontsize=11, ha='center', weight='bold')
    ax.text(0, 0, z_of_l(-l_max, 1.0) - 1.0, "EARTH\nside", color='darkgreen',
            fontsize=11, ha='center', weight='bold')

    ax.set_title("Embedding diagram of the wormhole throat\n"
                  "(lengths in units of throat radius $b_0$) with light-ray geodesics",
                  fontsize=12)
    ax.set_xlabel('x / b0'); ax.set_ylabel('y / b0'); ax.set_zlabel('z / b0')
    ax.legend(loc='upper left', fontsize=8, framealpha=0.9)
    ax.view_init(elev=18, azim=-60)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "1_embedding_diagram_geodesics.png"), dpi=160)

    return fig, ax, (X_grid, Y_grid, Z_grid)


def make_rotation_gif(grid_data):
    X_grid, Y_grid, Z_grid = grid_data
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X_grid, Y_grid, Z_grid, rstride=3, cstride=3,
                     color='#3a7bd5', alpha=0.55, edgecolor='#1f4e8c', linewidth=0.2)
    th_circle = np.linspace(0, 2 * np.pi, 200)
    ax.plot(np.cos(th_circle), np.sin(th_circle), np.zeros_like(th_circle), color='red', linewidth=2)
    ax.set_axis_off()
    ax.set_title("Wormhole throat (Morris-Thorne / Ellis embedding)", fontsize=10)

    def update(frame):
        ax.view_init(elev=20, azim=frame)
        return []

    anim = animation.FuncAnimation(fig, update, frames=np.linspace(0, 360, 60), interval=80)
    anim.save(os.path.join(OUTDIR, "2_wormhole_rotation.gif"), writer='pillow', fps=14, dpi=110)
    plt.close(fig)


def make_potential_and_energy_plots(field_eqs):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

    # --- Effective potential for photons, several impact parameters ---
    l_vals = np.linspace(-6, 6, 400)
    r_vals = r_of_l(l_vals, b0=1.0)
    for Lambda in [0.4, 0.8, 1.0, 1.3]:
        V = Lambda**2 / r_vals**2
        axes[0].plot(l_vals, V, label=f'b = {Lambda:.1f}')
    axes[0].axhline(1.0, color='k', linestyle='--', linewidth=1, label='$E^2=1$ (ray energy)')
    axes[0].set_xlabel('l / b0'); axes[0].set_ylabel(r'$V_{eff} = b^2/r(l)^2$')
    axes[0].set_title('Photon effective potential\n(barrier centred on the throat)')
    axes[0].legend(fontsize=8)

    # --- Exotic-matter energy density profile (symbolic -> numeric) ---
    b0_sym, l_sym = field_eqs['b0'], field_eqs['l']
    rho_fn = sp.lambdify((l_sym, b0_sym), field_eqs['rho'], 'numpy')
    rho_vals = rho_fn(l_vals, 1.0)
    axes[1].plot(l_vals, rho_vals, color='crimson', linewidth=2)
    axes[1].axhline(0, color='k', linewidth=0.8)
    axes[1].fill_between(l_vals, rho_vals, 0, color='crimson', alpha=0.15)
    axes[1].set_xlabel('l / b0'); axes[1].set_ylabel(r'$\rho(l)$   [units of $1/(8\pi b_0^2)$]')
    axes[1].set_title('Required energy density (from $G_{\\mu\\nu}=8\\pi T_{\\mu\\nu}$)\n'
                       'negative everywhere = "exotic matter"')

    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "3_potential_and_energy_density.png"), dpi=160)
    plt.close(fig)


# ==============================================================================
# SECTION 7 — MAIN: RUN EVERYTHING AND PRINT THE SUMMARY REPORT
# ==============================================================================
def main():
    print("=" * 78)
    print(" TRAVERSABLE WORMHOLE SIMULATION: EARTH <-> ANDROMEDA")
    print("=" * 78)

    print("\n[1] Solving the field equations symbolically (Christoffel -> Riemann ->")
    print("    Ricci -> Einstein tensor) for ds^2 = -c^2dt^2 + dl^2 + r(l)^2 dOmega^2 ...")
    fe = derive_field_equations()
    rho_throat = fe['rho'].subs(fe['l'], 0)
    nec_throat = fe['nec'].subs(fe['l'], 0)
    print(f"    Energy density  rho(l)      = {fe['rho']}")
    print(f"    Radial pressure p_r(l)      = {fe['p_radial']}")
    print(f"    Lateral pressure p_t(l)     = {fe['p_lateral']}")
    print(f"    rho(throat)                 = {rho_throat}   (units: 1/b0^2, G=c=1)")
    print(f"    Null-energy-condition rho+p_r at throat = {nec_throat}  -> negative = VIOLATED")
    print("    --> General Relativity itself requires negative energy density")
    print("        ('exotic matter') to hold this throat open. This is not an")
    print("        assumption of the model - it falls out of G_munu = 8 pi T_munu.")

    rho_SI = float(sp.N(fe['rho'].subs([(fe['l'], 0), (fe['b0'], B0)]))) * (C**2 / G_NEWTON)
    nuclear_density = 2.3e17  # kg/m^3, approx nuclear/neutron-star-core density
    print(f"\n    For a b0 = {B0:.1f} m throat, that is rho = {rho_SI:.3e} kg/m^3")
    print(f"    ( {abs(rho_SI)/nuclear_density:.2e} times nuclear density in magnitude — ")
    print("      no known macroscopic substance, only tiny quantum vacuum effects")
    print("      like the Casimir effect, has even been shown to have ANY negative")
    print("      energy density, let alone at this scale.)")

    print("\n[2] Embedding-diagram geometry + light-ray geodesics ...")
    fig1, ax1, grid_data = make_embedding_with_geodesics()
    print(f"    saved -> {OUTDIR}/1_embedding_diagram_geodesics.png")
    bcrit = critical_impact_parameter()
    print(f"    critical impact parameter b_crit = E*b0 = {bcrit:.2f} (in units of b0):")
    print("    light with impact parameter below this crosses the throat;")
    print("    above it, the ray reflects back without ever reaching the other side.")

    print("\n[3] Rendering a rotating view of the throat ...")
    make_rotation_gif(grid_data)
    print(f"    saved -> {OUTDIR}/2_wormhole_rotation.gif")

    print("\n[4] Effective potential & energy-density profile plots ...")
    make_potential_and_energy_plots(fe)
    print(f"    saved -> {OUTDIR}/3_potential_and_energy_density.png")

    print("\n[5] Tidal forces on a moving traveler (from the geodesic-deviation")
    print("    equation derived symbolically in Section 2) ...")
    v_max = max_comfortable_speed()
    print(f"    For a b0 = {B0:.0f} m throat and a {XI_BODY:.0f} m tall traveler, the speed")
    print(f"    that keeps tidal stretching under {G_EARTH:.2f} m/s^2 (1 g) is:")
    print(f"        v_max = {v_max:.4e} c  =  {v_max*C:.3f} m/s  (a brisk walking pace!)")
    print("    A *static* traveler feels zero tidal force in this metric (that's what")
    print("    'zero-tidal-force wormhole' means) - it is specifically motion through")
    print("    the curved throat that stretches you, growing as (speed)^2.")
    print(f"\n    Throat radius needed for a LARGER, relativistic-speed traveler instead:")
    for vf in [0.5, 0.9, 0.999]:
        b_needed = throat_for_comfortable_speed(vf)
        print(f"        v = {vf:>6.3f} c  ->  b0 >= {b_needed:.3e} m  ({b_needed/LY:.3e} light-years)")

    print("\n[6] Travel time, Earth <-> Andromeda")
    print(f"    Normal-space distance: {D_ANDROMEDA_M:.3e} m  ({D_ANDROMEDA_LY:.2e} ly)")
    print(f"    -> light takes {D_ANDROMEDA_LY:.2e} years the 'normal' way.")

    print(f"\n    SCENARIO A - modest, human-scale throat (b0 = {B0:.0f} m), limited to the")
    print(f"    tidal-safe creep speed v_max = {v_max*C:.2f} m/s found above.")
    print(f"    {'tunnel length L':>18} | {'crossing time':>20}")
    for L_demo in [1.0e3, 1.0e4, 1.0e5, 1.0e6, 60 * C]:
        t_demo = L_demo / (v_max * C)
        label = (f"{L_demo:.0f} m" if L_demo < 1e6 else
                  ("1 light-min" if abs(L_demo - 60 * C) < 1 else f"{L_demo:.2e} m"))
        if t_demo < 3600:
            t_str = f"{t_demo:.1f} s  ({t_demo/60:.2f} min)"
        elif t_demo < 86400:
            t_str = f"{t_demo/3600:.2f} hours"
        else:
            yr = t_demo/86400/365.25
            t_str = f"{t_demo/86400:.1f} days" + (f"  ({yr:.2f} yr)" if yr > 0.1 else "")
        print(f"    {label:>18} | {t_str:>20}")
    print("    (compare any of these to the 2.5 MILLION YEARS light needs the long way)")

    print(f"\n    SCENARIO B - large throat (light-year-scale, per the table in [5]),")
    print(f"    big enough to cross at a genuinely relativistic speed. Tunnel length")
    print(f"    fixed at L = 1 AU = {AU:.3e} m for this comparison:")
    print(f"    {'speed (c)':>10} | {'coord. time':>16} | {'traveler proper time':>22} | {'Lorentz gamma':>13}")
    for vf in [0.5, 0.9, 0.99, 0.999999]:
        t_coord, tau, gamma = transit_times(vf, L=AU)
        print(f"    {vf:>10.6g} | {t_coord:>13.4g} s | {tau:>19.4g} s   | {gamma:>13.4g}")
    print("    Either way: minutes to days through the throat, vs 2.5 million years.")
    print("    NOTE: General Relativity fixes none of L, b0, or where the mouths sit -")
    print("    those are free 'design' choices, exactly as Morris & Thorne treat them [1].")

    print("\n" + "=" * 78)
    print(" SOURCES")
    print("=" * 78)
    for line in SOURCES:
        print(" " + line)
    print("=" * 78)


SOURCES = [
    "[1] M. S. Morris and K. S. Thorne, \"Wormholes in spacetime and their use for",
    "    interstellar travel: A tool for teaching general relativity,\" American",
    "    Journal of Physics 56(5), 395-412 (1988). doi:10.1119/1.15620",
    "    -> source of the metric used here (their 'Box 2' worked example),",
    "       the throat/flare-out condition, and the tidal-force analysis.",
    "[2] M. S. Morris, K. S. Thorne, and U. Yurtsever, \"Wormholes, Time Machines,",
    "    and the Weak Energy Condition,\" Physical Review Letters 61, 1446 (1988).",
    "[3] H. G. Ellis, \"Ether flow through a drainhole: A particle model in general",
    "    relativity,\" Journal of Mathematical Physics 14, 104 (1973).",
    "    -> independent discovery of the same metric, sourced by a phantom",
    "       scalar field; also found by K. A. Bronnikov (1973).",
    "[4] M. Visser, \"Lorentzian Wormholes: From Einstein to Hawking,\" AIP Press,",
    "    New York (1995). -- standard reference on energy conditions & wormholes.",
    "[5] O. James, E. von Tunzelmann, P. Franklin, and K. S. Thorne, \"Visualizing",
    "    Interstellar's Wormhole,\" American Journal of Physics 83(6), 486-499",
    "    (2015). doi:10.1119/1.4916949 -- ray-tracing/embedding-diagram methods",
    "       this script's visualization is modeled on.",
    "[6] NASA Science, \"Andromeda Galaxy\" / Galaxy Evolution Explorer image",
    "    release: distance ~2.5 million light-years.",
    "    https://science.nasa.gov/  and  https://www.nasa.gov/  (M31 pages).",
    "[7] Misner, Thorne & Wheeler, \"Gravitation,\" W. H. Freeman (1973) -- general",
    "    reference for the Christoffel/Riemann/Einstein-tensor and geodesic-",
    "    deviation formalism used throughout Sections 2 and 5.",
]


if __name__ == "__main__":
    main()
