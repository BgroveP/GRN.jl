# Third-party notices

GRN.jl is licensed under the MIT License (see `LICENSE`). Parts of it are
adapted from the third-party projects listed below, which are distributed
under their own licenses. Their copyright notices and license texts are
reproduced here as those licenses require.

---

## posterior (R package)

- **URL:** <https://github.com/stan-dev/posterior>
- **License:** BSD 3-Clause
- **Used in:** `src/diagnostics.jl` — `_ess`, `_split_chains`,
  `_rank_normalize`, and the `kind=:bulk` / `:tail` / `:basic` logic of `ess`
  are adapted from the effective sample size functions of `posterior`
  (`ess_bulk`, `ess_tail`, `ess_basic`, and the internal `.ess`), which
  implement Vehtari et al. (2021), *Bayesian Analysis* 16(2), 667–718.
  Deviations: ranks and tail quantiles use all draws before the chains are
  split (as in MCMCDiagnosticTools.jl), and fewer than 5 draws per split
  chain return `NaN`.

```
Copyright (c) 2021, posterior package authors; Stan Developers and their
Assignees; Trustees of Columbia University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

---

## Published algorithms (no license required)

Implemented from the publications; listed for reference.

- `_qnorm` (`src/diagnostics.jl`): Wichura, M. J. (1988). Algorithm AS 241:
  The percentage points of the normal distribution. *Journal of the Royal
  Statistical Society, Series C (Applied Statistics)* 37(3), 477–484.

## Used for validation only (no code copied)

- MCMCDiagnosticTools.jl (MIT), <https://github.com/TuringLang/MCMCDiagnosticTools.jl>:
  reference values in the "Effective sample size" testset of
  `test/runtests.jl`.
