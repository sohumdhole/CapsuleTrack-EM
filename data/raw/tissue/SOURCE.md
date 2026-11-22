# Tissue dielectric parameters — provenance

`gabriel1996_colecole.csv` holds the four-term Cole–Cole fitting parameters of
Gabriel, Lau & Gabriel (1996) for the abdominal and reference tissues used in this
project.

* Primary reference: S. Gabriel, R. W. Lau, C. Gabriel, "The dielectric properties of
  biological tissues: III. Parametric models for the dielectric spectrum of tissues",
  *Phys. Med. Biol.* 41 (1996) 2271–2293. doi:10.1088/0031-9155/41/11/003
* Transcribed from the IFAC-CNR public copy of the Gabriel report, Appendix C
  ("Summary of fitting parameters"):
  https://niremf.ifac.cnr.it/docs/DIELECTRIC/AppendixC.html
* Units: tau1 in ps, tau2 in ns, tau3 in µs, tau4 in ms, sig (ionic conductivity) in S/m.
* The IFAC page prints conductivity to three decimals; the skin ionic conductivities
  (0.0002 / 0.0004 S/m) are the values of the original report.

Model (Gabriel et al. 1996, eq. 3):

    eps*(w) = ef + sum_n  del_n / (1 + (j w tau_n)^(1 - alf_n))  +  sig / (j w eps0)

The IT'IS Foundation Tissue Properties Database (V4.x) builds its low-frequency
conductivity table on these same Gabriel dispersion relations plus later literature.
The IT'IS download host was not reachable from the build environment, so the Gabriel
parameters are used directly. Swapping in the IT'IS V4.2 LF-conductivity table only
requires editing `capsuletrack/tissue.py::lf_conductivity_override`.

Validation check (tests/test_tissue.py): muscle conductivity at 10 kHz from this table
evaluates to ~0.34 S/m, consistent with the published Gabriel curves.
