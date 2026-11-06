# Third-party code: Failenschmid et al. (2025)

The files classes.R and functions.R in this folder are copied unmodified
from the supplementary material of:

Failenschmid, J. I., Vogelsmeier, L. V. D. E., Mulder, J., & Jongerling, J.
(2025). Modeling non-linear psychological processes: Reviewing and evaluating
non-parametric approaches and their applicability to intensive longitudinal
data. British Journal of Mathematical and Statistical Psychology.
https://doi.org/10.1111/bmsp.12397

Copyright (c) 2024 Jan Ian Failenschmid. Licensed under the GNU General
Public License v3.0 or later (see gpl-3.0.txt).

They are sourced by ../power_sim.R, which uses the gen_model class
and the functions add_conditions(), sim_tsm() and get_tsm_data() to simulate
data from the generative models (exponential growth, logistic growth, damped
oscillator, cusp catastrophe). The simulate() function and the model
definitions in power_sim.R are adapted from the same source.