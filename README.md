# Figure Code for the Midterm Project

This folder contains the MATLAB code needed to regenerate the figures used in the final manuscript and appendix.

## Files

- `plot_figures_B.m`: shared generator. Run it with no input to regenerate all figures, or with a figure name to regenerate one figure.
- `Fig1_boundary_sampling.m`: entry script for Fig. 1.
- `Fig2_finite_size_diagnostics.m`: entry script for Fig. 2.
- `Fig3_tamm_bound_state.m`: entry script for Fig. 3.
- `Fig4_ssh_bands_edge.m`: entry script for Fig. 4.
- `FigS1_ldos_map.m`: entry script for Appendix Fig. S1.
- `FigS2_ipr_splitting.m`: entry script for Appendix Fig. S2.
- `compute_finite_size_scaling.m`: helper for the log-log finite-size IDOS scaling panel in Fig. 2.

## How to Run

From this folder, run:

```matlab
plot_figures_B
```

To regenerate a single figure, run the matching script, for example:

```matlab
Fig2_finite_size_diagnostics
FigS1_ldos_map
```

The output files will be written to:

```text
figures_B/
```

## Figures

- `Fig1_boundary_sampling`: electronic and phonon boundary quantization sampling.
- `Fig2_finite_size_diagnostics`: IDOS scaling, broadened spectra, and separated electronic/phonon spacing distributions.
- `Fig3_tamm_bound_state`: Tamm-like surface bound state and exponential localization.
- `Fig4_ssh_bands_edge`: SSH bulk gap and open-chain in-gap edge states.
- `FigS1_ldos_map`: optional LDOS localization diagnostic.
- `FigS2_ipr_splitting`: optional IPR and finite-size edge-state splitting diagnostic.

The final manuscript uses Fig. 1--Fig. 4 in the main text, so the main text stays below the five-figure limit. Fig. S1--Fig. S2 are appendix diagnostics that document the localization workflow.
