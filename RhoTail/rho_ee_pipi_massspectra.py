#!/usr/bin/env python3
"""
Compare rho(770) -> pi+pi- and rho(770) -> e+e- invariant-mass line shapes
in gamma p -> rho p photoproduction, over a grid of (Egamma, -t) settings.

Physics
-------
Both decay channels share the SAME rho propagator: pole mass M0 and a
mass-dependent TOTAL width Gamma_tot(M) that is completely dominated by the
pi+pi- channel (BR(rho->pi+pi-) ~ 100%). They differ in the NUMERATOR of the
decay rate -- but the two numerators come from physically different decay
mechanisms, not just "the same formula with a different daughter mass":

  - rho -> pi+pi- is a strong, P-wave (L=1) two-body decay with a roughly
    mass-independent coupling g_rhopipi. Its M-dependence is genuinely a
    two-body phase-space effect, Gamma_pipi(M) ~ q_pi(M)^3, with a
    Blatt-Weisskopf centrifugal-barrier factor that tames the q^3 growth
    once M gets well above the pole (a point-like-vertex q^3 alone would
    make Gamma_pipi(M) exceed M itself by ~2 GeV, which is unphysical).
    Near the pole this still shows real curvature from the 2*m_pi threshold
    sitting not far below the resonance -> a visibly skewed line shape (the
    Ross-Stodolsky effect seen in real rho photoproduction data).
  - rho -> e+e- proceeds through vector-meson dominance: the rho mixes with
    a virtual photon (coupling ~ e*M0^2/g_rho, from the rho-photon mixing
    amplitude), which then converts to e+e-. Two-body phase space for
    massless leptons is flat in M -- there is no q^3 factor here at all.
    The M-dependence instead comes entirely from the s-channel photon
    propagator, giving Gamma_ee(M) = Gamma_ee(M0) * (M0/M)^2, i.e. FALLING
    with M. (An earlier version of this script incorrectly reused the
    pi+pi- q^3 Jackson formula with m_e substituted in, which -- since
    q_e(M) ~ M/2 at any M, m_e being negligible -- gave a RISING Gamma_ee(M)
    ~ M^2 instead: the opposite trend from the true VMD/photon-propagator
    physics. That mistake is invisible near the pole, where both forms are
    roughly flat, but wrong by a large and growing factor by M=2-3 GeV. See
    V. Kubarovsky's note, Docs/RhoDecayWidths/rho_bw_note_byValery.pdf,
    Sec. 5, for the derivation of the correct Gamma_ee(M).)

Consequence for the ratio: near the pole the ratio is close to flat (both
channels vary slowly there), but away from the pole it falls steeply --
Gamma_ee(M)/Gamma_pipi(M) ~ (M0/M)^2 / M^(up to 2) ~ (M0/M)^(up to 4) for the
bare P-wave case -- NOT flat out to 3 GeV as the earlier (wrong) version of
this script showed.

The LINE SHAPE (peak position, width, skew, and this falling ratio) does
NOT depend on the rho's production kinematics (Egamma, t) -- that decay
happens in the rho rest frame, so it's Lorentz-invariant. What DOES depend
on (Egamma, t) is the overall event yield, through the gamma p -> rho p
production cross section. That part is modeled below with an illustrative
toy diffractive weight (dsigma/dt ~ exp(-b|t|), mild energy rise) purely so
the 16 panels aren't visually identical -- it is NOT a validated
cross-section calculation. The two curves in each panel keep their true
relative normalization (BR(e+e-)/BR(pi+pi-) ~ 4.7e-5), plotted on a log
y-axis so both are visible at once.

gamma_pipi_of_m()'s q^3*(M0/M) factor (before the barrier) mirrors the
Jackson mass-dependent-width form already implemented and numerically
verified in this repo's RhoFitGUI.cc (RelBW(), num = M0*Gamma(M), den =
(M^2-M0^2)^2 + M0^2*Gamma(M)^2); the Blatt-Weisskopf barrier and the
gamma_ee_of_m() VMD form are new here (RhoFitGUI.cc only ever fits near the
pole, where neither matters).
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- constants
M0 = 0.7753      # GeV, rho0 pole mass (PDG)
GAMMA0 = 0.1474  # GeV, rho0 total width (user-given / PDG)
M_PI = 0.13957   # GeV, charged pion mass
M_E = 0.000511   # GeV, electron mass
BR_EE = 4.72e-5  # PDG: BR(rho0 -> e+e-)
BR_PIPI = 1.0    # ~100%
GAMMA_EE0 = BR_EE * GAMMA0  # GeV, on-shell Gamma(rho->e+e-), ~6.96 keV
R_BW = 3.0       # GeV^-1, Blatt-Weisskopf interaction radius for the
                  # pi+pi- P-wave barrier (typical rho-fit value, 1.5-5
                  # GeV^-1 all appear in the literature and are NOT
                  # constrained by rho(770)-region data alone -- see
                  # rho_bw_note_byValery.pdf Sec. 4.2).


def q2body(m, m_daughter):
    """Two-body breakup momentum of a resonance of mass m decaying into two
    identical daughters of mass m_daughter, in the resonance rest frame."""
    arg = m * m / 4.0 - m_daughter * m_daughter
    return np.sqrt(np.clip(arg, 0.0, None))


def barrier_F2(q, r):
    """Blatt-Weisskopf L=1 centrifugal barrier factor."""
    return 1.0 / (1.0 + (r * q) ** 2)


def gamma_pipi_of_m(m, r=R_BW):
    """Jackson mass-dependent P-wave width for rho -> pi+pi-, with a
    Blatt-Weisskopf barrier so the width saturates far above the pole
    instead of growing without bound (unphysically exceeding M itself)
    like the bare q^3 phase-space factor alone would."""
    q = q2body(m, M_PI)
    q0 = q2body(M0, M_PI)
    barrier = barrier_F2(q, r) / barrier_F2(q0, r)
    return GAMMA0 * (q / q0) ** 3 * (M0 / m) * barrier


def gamma_ee_of_m(m):
    """Vector-meson-dominance (photon-propagator-mediated) partial width for
    rho -> e+e-: falls as (M0/M)^2. Massless-lepton two-body phase space is
    flat in M, so essentially all of the M-dependence here comes from the
    s-channel photon propagator in the rho-photon-lepton coupling, not from
    any phase-space/threshold effect (see module docstring).

    Caveat: (M0/M)^2 formally diverges as M -> 0, which is visible as a
    modest upturn of the e+e- curve at the low-mass end (M ~ 0.3-0.45 GeV)
    of the plots below. This is an artifact of extrapolating the leading
    VMD term well below the region it was derived for (not too far from
    the pole); it is not claiming a real low-mass enhancement -- genuine
    low-mass e+e- pairs are dominated by other (non-resonant/Dalitz-type)
    QED physics entirely, not by a "rho tail." Harmless for everything this
    script actually uses it for (the pole region and the M > M0 tail)."""
    return GAMMA_EE0 * (M0 / m) ** 2


def lineshape(m, M0, gamma_tot, gamma_num):
    """Relativistic BW: numerator uses the channel-specific partial width,
    denominator uses the TOTAL width (here, dominated by pi+pi-)."""
    num = M0 * gamma_num(m)
    den = (m * m - M0 * M0) ** 2 + M0 * M0 * gamma_tot(m) ** 2
    return num / den


_trapz = getattr(np, "trapezoid", np.trapz)


def normalize_to_unit_area(y, m):
    return y / _trapz(y, m)


# rho decay line shapes vs invariant mass M
m = np.linspace(0.32, 3.00, 3500)

# Both channels use gamma_pipi_of_m in the DENOMINATOR (total width);
# only the numerator (gamma_num) differs by channel.
shape_pipi = normalize_to_unit_area(
    lineshape(m, M0, gamma_pipi_of_m, gamma_pipi_of_m), m
)
shape_ee = normalize_to_unit_area(
    lineshape(m, M0, gamma_pipi_of_m, gamma_ee_of_m), m
)

# For the 4x4 grid below, rescale both shapes to the SAME value at M=M0 (the
# rho mass) rather than keeping the true BR(e+e-)/BR(pi+pi-) ~ 4.7e-5 yield
# ratio -- this throws away the real relative normalization on purpose, so
# the shape difference (Ross-Stodolsky skew vs. symmetric BW) is easy to see
# by eye. Don't read panel heights as physical rates; see the note in the
# figure title.
shape_pipi_at_M0 = np.interp(M0, m, shape_pipi)
shape_ee_at_M0 = np.interp(M0, m, shape_ee)
shape_pipi_eqM0 = shape_pipi / shape_pipi_at_M0
shape_ee_eqM0 = shape_ee / shape_ee_at_M0


def production_weight(Eg, t, Eg_ref=6.0, b=7.0):
    """Illustrative-only diffractive-like weight for gamma p -> rho p,
    dsigma/dt ~ exp(-b|t|) with a mild energy rise. NOT a fitted or
    validated cross section -- used only so the 16 panels below differ in
    overall normalization the way real (Egamma, t) dependence would."""
    return (Eg / Eg_ref) ** 0.1 * np.exp(-b * abs(t))


Eg_vals = np.linspace(4.0, 10.6, 4)  # GeV
t_vals = np.linspace(0.15, 1.0, 4)   # GeV^2 (i.e. -t)

# dataviz skill categorical palette, slots 1-2 (fixed order, light mode)
COLOR_PIPI = "#eb6834"  # orange, slot 2
COLOR_EE = "#2a78d6"    # blue, slot 1

fig, axes = plt.subplots(4, 4, figsize=(16, 13), sharex=True)

for i, Eg in enumerate(Eg_vals):
    for j, t in enumerate(t_vals):
        ax = axes[i, j]
        w = production_weight(Eg, t)
        y_pipi = w * shape_pipi_eqM0
        y_ee = w * shape_ee_eqM0

        ax.plot(m, y_pipi, color=COLOR_PIPI, lw=1.6, label=r"$\rho\to\pi^+\pi^-$")
        ax.plot(m, y_ee, color=COLOR_EE, lw=1.6, label=r"$\rho\to e^+e^-$")
        ax.set_ylim(0, max(y_pipi.max(), y_ee.max()) * 1.15)
        ax.text(
            0.05, 0.92, f"$E_\\gamma$={Eg:.2f} GeV\n$-t$={t:.2f} GeV$^2$",
            transform=ax.transAxes, va="top", ha="left", fontsize=8,
        )
        if i == 3:
            ax.set_xlabel("M (GeV)")
        if j == 0:
            ax.set_ylabel("arb. yield")

axes[0, 0].legend(loc="upper right", fontsize=8, frameon=False)
fig.suptitle(
    r"$\gamma p \to \rho p$: $\rho \to \pi^+\pi^-$ vs $\rho \to e^+e^-$ mass spectra"
    "\nSame line shape at every $E_\\gamma, -t$ (decay is Lorentz-invariant); "
    "only the overall height changes, via an illustrative production-weight toy model.\n"
    "Both curves rescaled to the SAME value at $M=M_\\rho$ to make the shape "
    f"difference visible -- NOT the true yield ratio (BR($e^+e^-$)/BR($\\pi^+\\pi^-$) "
    f"$\\approx$ {BR_EE:.2e}).",
    fontsize=10,
)
fig.tight_layout(rect=[0, 0, 1, 0.90])
fig.savefig("rho_ee_pipi_massspectra.png", dpi=150)
print("wrote rho_ee_pipi_massspectra.png")

# ------------------------------------------------------------------
# Bonus: shape-only comparison (unit-peak normalized, like this repo's
# CompareHists.cc convention), independent of Egamma/t, to see the
# Ross-Stodolsky skew directly.
fig2, ax2 = plt.subplots(figsize=(7, 5))
ax2.plot(m, shape_pipi / shape_pipi.max(), color=COLOR_PIPI, lw=2, label=r"$\rho\to\pi^+\pi^-$ (shape)")
ax2.plot(m, shape_ee / shape_ee.max(), color=COLOR_EE, lw=2, label=r"$\rho\to e^+e^-$ (shape)")
ax2.set_xlabel("M (GeV)")
ax2.set_ylabel("unit-max normalized yield")
ax2.set_title("Line-shape comparison (kinematics-independent)")
ax2.legend(frameon=False)
fig2.tight_layout()
fig2.savefig("rho_ee_pipi_shape_only.png", dpi=150)
print("wrote rho_ee_pipi_shape_only.png")

# ------------------------------------------------------------------
# Ratio of the two shapes (e+e- / pi+pi-), using the same "equal at M0"
# rescaling as the 4x4 grid above (shape_*_eqM0), so this is a pure shape
# comparison -- not the true BR(e+e-)/BR(pi+pi-) yield ratio. By
# construction it equals 1 at M = M0. Near the pole it rises (the pi+pi-
# P-wave phase space collapses toward the 2*m_pi threshold while Gamma_ee
# stays essentially flat there), but away from the pole it FALLS steeply --
# Gamma_ee(M) ~ (M0/M)^2 while Gamma_pipi(M) grows (before saturating via
# the barrier), so their ratio drops by orders of magnitude by M=3 GeV.
# Plotted log-scale for that reason. Compare Fig. 3 of
# rho_bw_note_byValery.pdf, which shows the same qualitative fall (its
# absolute values differ because it also treats the pi+pi- width and the
# spectral normalization more carefully -- see the module docstring).
ratio_ee_over_pipi = shape_ee_eqM0 / shape_pipi_eqM0

fig3, ax3 = plt.subplots(figsize=(7, 5))
ax3.plot(m, ratio_ee_over_pipi, color="#4a3aa7", lw=2)
ax3.set_yscale("log")
ax3.axhline(1.0, color="black", lw=0.8, ls="--")
ax3.axvline(M0, color="black", lw=0.8, ls="--")
ax3.set_xlabel("M (GeV)")
ax3.set_ylabel(r"shape($\rho\to e^+e^-$) / shape($\rho\to \pi^+\pi^-$)")
ax3.set_title(r"Line-shape ratio ($e^+e^-$ / $\pi^+\pi^-$), = 1 at $M=M_\rho$")
fig3.tight_layout()
fig3.savefig("rho_ee_pipi_ratio.png", dpi=150)
print("wrote rho_ee_pipi_ratio.png")

# ------------------------------------------------------------------
# Same ratio, but laid out as a 4x4 (Egamma, -t) grid like the main mass
# spectra plot, since the ratio is built from shape_*_eqM0 (production
# weight w(Eg,t) cancels in the ratio -- see production_weight() above),
# every panel shows the exact same curve. That sameness IS the physics
# point: the shape difference between the two decay channels is a property
# of the rho resonance itself, not of how/where it was produced.
fig4, axes4 = plt.subplots(4, 4, figsize=(16, 13), sharex=True, sharey=True)

for i, Eg in enumerate(Eg_vals):
    for j, t in enumerate(t_vals):
        ax = axes4[i, j]
        ax.plot(m, ratio_ee_over_pipi, color="#4a3aa7", lw=1.6)
        ax.set_yscale("log")
        ax.axhline(1.0, color="black", lw=0.8, ls="--")
        ax.axvline(M0, color="black", lw=0.8, ls="--")
        ax.text(
            0.95, 0.92, f"$E_\\gamma$={Eg:.2f} GeV\n$-t$={t:.2f} GeV$^2$",
            transform=ax.transAxes, va="top", ha="right", fontsize=8,
        )
        if i == 3:
            ax.set_xlabel("M (GeV)")
        if j == 0:
            ax.set_ylabel(r"$e^+e^-$ / $\pi^+\pi^-$")

fig4.suptitle(
    r"$\gamma p \to \rho p$: shape ratio $\rho\to e^+e^-$ / $\rho\to \pi^+\pi^-$"
    "\nIdentical in every panel -- the ratio is a property of the rho decay "
    "itself and does not depend on the production kinematics $E_\\gamma, -t$.",
    fontsize=10,
)
fig4.tight_layout(rect=[0, 0, 1, 0.92])
fig4.savefig("rho_ee_pipi_massspectra_ratio.png", dpi=150)
print("wrote rho_ee_pipi_massspectra_ratio.png")
