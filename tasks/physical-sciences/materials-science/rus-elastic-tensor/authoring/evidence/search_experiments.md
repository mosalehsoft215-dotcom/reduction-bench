# Global search: what works on these sweeps

All runs used the 59 merged resonances produced by `solution/rus_invert.py`, the environment image, and one CPU per run. "Right" means every sample-frame constant within 3% of the true one.

| strategy | right basin | notes |
|---|---|---|
| 196 Sobol starts, each a damped Gauss-Newton with Hungarian re-matching (order 8, lowest 20, caps 6% down to 1.5%) | 0 / 196 | local fits stall in whatever basin they start in |
| DE, linear bounds, order 8, lowest 20, cap 2% (popsize 24, 180 generations) | 2 / 8 | wrong basins score 0.41-0.55 against 0.33 for the right one, too close to trust |
| DE, log-space diagonals, C12'/C13' as fractions of the Born limit, otherwise as above | 2 / 8 | |
| DE, same coordinates, order 10, lowest 30, cap 1% | 1 / 8 | clean separation (right 1.33, wrong at least 9.0) but a much rougher landscape |
| continuation DE: (order 8, 12 res, 4%) then (8, 20, 2%) then (10, 30, 1%), each stage seeded with the last population | 4 / 8 | right 0.351, wrong at least 8.93; this is what the reference solution uses |

With a 50% hit rate per seed and an unambiguous score, the reference runs seeds two at a time and stops at the first clear success. Twelve seeds without one success has odds of about 1 in 4000. On this machine seed 1 succeeds in the first round.

The scripts for these runs are in `search/`. They expect the bundle mounted at `/b` inside the environment image, and they take the strategy letter and a comma-separated seed list as arguments:

    docker run --rm -v <bundle>:/b -v <this dir>/search:/s <env-image> python /s/de_variants.py G 1,2,3,4,5,6,7,8
