# Physical scale and applicability of the synthetic case

The design draws on published pulsed-electrosynthesis workflows described in
`sources.md`. Its concentrations, kinetic constants, blocking law and operating
scenarios are declared case assumptions. No particular molecule, measured rate
constant, commercial cell or completed laboratory experiment is represented.

The five liquid layers are 85-110 micrometres thick. The nominal layer is 100
micrometres, the substrate concentration is 5 mmol/L, and the nominal diffusion
coefficients are 1.0e-5 and 0.8e-5 cm^2/s. These inputs define a diffusion length
sqrt(Ds*T) of approximately 77 micrometres over the six-second horizon. The
dimensionless diffusion time Ds*T/L^2 is 0.6 in the nominal cell. Spatial profiles
therefore have time to evolve substantially without instantaneously equilibrating
through the layer. These are dimensional consequences of the case inputs, not
measurements from the cited papers.

The nominal substrate inventory is C0*L/1000 = 5e-8 mol/cm^2. The public charge cap
of 0.70 electron equivalents corresponds to approximately 3.38 mC/cm^2, or a mean
cathodic current density of 0.563 mA/cm^2 over six seconds. A cathode area is not
required for grading because all inventories and charge are normalized consistently
per unit area. Converting them to an actual device current would require its area.

At R=C0 the initial homogeneous consumption timescale is
1/(2*kd*C0+kh), approximately 0.048 s in the nominal case. It lengthens at lower
intermediate concentrations. Nominal surface-coverage relaxation is about 0.17 s
at -1.30 V and 0.10 s at -0.80 V. These scales make intermediate accumulation,
reaction and electrode recovery compete over pulses of tens to hundreds of
milliseconds. They explain why a time-dependent scientific model is needed.

The model intentionally omits charging, uncompensated resistance, migration,
convection, molecular side reactions beyond the declared network, and spatial
transport of inert products. The input voltage is the ideal interfacial potential.
Applying a computed program experimentally would require validation of those
approximations, identification of the chemistry and its rates, and conversion to
the actual potentiostat's method format. The supplied waveform limits are case
design requirements, not verified specifications of a named instrument.

The computational task is fully determined by local equations and inputs. Its
scientific scope is robust experiment design for a disclosed reduced model; experimental
predictiveness and a global optimum are not asserted.

For the nominal inventory, the 0.00025-equivalent charge-bin noise corresponds
to about 1.21 microcoulomb/cm^2, or 24.1 microampere/cm^2 when averaged over a
50 ms bin. The independent 0.0002-equivalent baseline prior corresponds to
0.965 microcoulomb/cm^2 per bin. The gain prior is 2%, endpoint assay noise is
one percentage point of initial substrate equivalents, and physical nuisance
priors have log standard deviations of 0.20-0.30. These are declared experimental
planning assumptions. They are not measurements or vendor specifications.
The information calculation is local at each of the five nominal scenarios.
