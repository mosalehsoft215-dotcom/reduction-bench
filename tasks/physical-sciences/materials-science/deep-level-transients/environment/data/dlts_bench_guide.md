# DLTS bench guide

This bench characterises electron traps in n-type silicon by deep-level transient spectroscopy. Each diode is an Au
Schottky dot evaporated through a 1.00 mm mask onto an n-type Si epitaxial layer grown on an n+ Si substrate, with a
back contact on the substrate. One diode folder holds everything recorded on one diode.

## The bench

- **Cryostat.** Closed-cycle, sample clamped to a copper stage. `T_K` in every file is the stage sensor reading
  after the stage has settled at that set point.
- **Capacitance meter.** It applies a 30 mV test signal at the frequency in `diode.json` (1 MHz) and reports the
  admittance between its two terminals as a capacitance `C_pF` in parallel with a conductance `G_uS`.
- **Pulse generator.** It holds the diode at the reverse bias and, once per period, steps it to the pulse voltage for
  the pulse width. Rise and fall times are below 2 ns.
- **Digitiser.** It records the meter output after the end of each pulse on a logarithmic time base from 20 us to
  2 s (256 points; 128 in the capture series), averaging 64 transients. The meter needs about 50 us to recover from
  the pulse, so the first points carry its recovery.

## Files in a diode folder

- `diode.json`: dot area (`area_cm2`), meter frequency, and the settings of the three transient series.
- `log.txt`: the operator's notes.
- `cv_300k.csv`: `bias_v, C_pF, G_uS`, from 0 to -10 V in 0.1 V steps at 300 K, each point taken in steady state.
- `scan.csv`: `T_K, t_s, C_pF, G_uS`. One transient per set point from 100 to 340 K in 2.5 K steps, at the scan's
  reverse bias and pulse settings. `t_s` is the time after the end of the pulse.
- `bias.csv`: `T_K, reverse_bias_v, t_s, C_pF, G_uS`. At five set points, transients at reverse biases of -2, -5
  and -10 V with the same pulse voltage and width.
- `capture.csv`: `T_K, pulse_width_s, t_s, C_pF, G_uS`. At two set points, transients after pulses of 13 widths
  from 5 ns to 50 us.
- `checks.csv`: `T_K, VF_10uA_v, VF_100uA_v`. After the transient at each scan set point, the forward voltage at
  10 uA and at 100 uA, the bench's routine contact check.

Rows of one transient are consecutive and in increasing `t_s`.

## Device physics used on this bench

Constants: q = 1.602176634e-19 C, k = 8.617333262e-5 eV/K, eps_s = 11.7 eps_0 with eps_0 = 8.8541878128e-14 F/cm,
N_c(T) = 2.86e19 (T/300)^1.5 cm^-3, v_th(T) = 2.30e7 (T/300)^0.5 cm/s. Depth x is measured from the metal.

**Depletion layer.** The depletion approximation holds: the layer between the metal and depth W is free of
electrons, and the junction capacitance is C = eps_s A / W for a dot of area A. With N+(x) the positive space-charge
density in the layer, the bias V fixes W through

    (q / eps_s) * integral_0^W x N+(x) dx = V_bi - V - kT/q,     V_bi = phi_B - (E_c - E_F),

with phi_B the barrier height of the contact and E_c - E_F = kT ln(N_c / n) in the neutral epi, n its free-electron
density. N+ is the donor density N_d(x) corrected for the charge held by the traps.

**The trap.** One electron trap of density N_T at energy E_T below the conduction band. A trap is either
acceptor-like (neutral when empty, negative when it holds an electron) or donor-like (positive when empty, neutral
when it holds an electron). In the neutral epi every trap holds an electron.

**Lambda layer.** Near the edge of the depletion layer, the traps lie below the Fermi level over a distance
lambda = sqrt(2 eps_s (E_F - E_T) / (q^2 N_d)), with N_d taken at the edge (lambda = 0 if E_T lies above E_F).
These traps stay filled at any bias. At reverse bias, the traps that can empty are those between W_p - lambda_p and
W_r - lambda_r, where W_p is the depletion width at the pulse voltage, W_r the steady width at the reverse bias, and
each lambda is evaluated at its own edge. Traps shallower than W_p - lambda_p are never refilled by the pulse.

**Rates.** A filled trap emits its electron at the rate

    e_n = sigma_n v_th(T) N_c(T) exp(-(E_c - E_T) / kT),

and an empty trap in neutral material captures one at c_n = sigma_n v_th(T) n, where sigma_n is the electron
capture cross-section. During a pulse of width t_p the traps in the window fill as
1 - exp(-c_n t_p); between pulses they empty as exp(-e_n t). The pulses repeat at the period in `diode.json`.

**Forward current.** Over the barrier, I = A A* T^2 exp(-q phi_B / kT) [exp(q V_j / (n kT)) - 1] with
A* = 112 A cm^-2 K^-2 for n-type Si, V_j the voltage across the junction and n the ideality factor of the contact.

## The values reported for a diode

- `level_ev`: E_c - E_T of the trap in zero electric field, in eV.
- `density_cm3`: N_T, in cm^-3.
- `sigma_300k_cm2`: the electron capture cross-section at 300 K, in cm^2.
