# fluctuation-analysis-skill

Public mirror of the **`fluctuation-analysis`** skill — turn Luria–Delbrück fluctuation-assay colony counts into mutation rates, compare rates, size experiments, and optimize a mutator induction window with hard protection caps.

- **Version:** 0.7.5
- **License:** MIT (see `LICENSE` and `fluctuation-analysis/LICENSE`)
- **Language:** Chinese docs; commands are language-neutral

## Where is the skill?

The skill itself lives in **[`fluctuation-analysis/`](fluctuation-analysis/)**.

This matters: DSH / Agent Skills loaders require the skill **directory name to match the `name` in `SKILL.md` frontmatter**. Keeping the skill in a subdirectory means a plain `git clone` gives you a correctly named bundle.

## Install

```bash
git clone https://github.com/onetree3636/fluctuation-analysis-skill.git
cp -r fluctuation-analysis-skill/fluctuation-analysis ~/.agents/skills/
# or symlink:
ln -s "$PWD/fluctuation-analysis-skill/fluctuation-analysis" ~/.agents/skills/fluctuation-analysis
```

For Claude Code, use `~/.claude/skills/` instead of `~/.agents/skills/`.

## Quick start

```bash
cd fluctuation-analysis
python3 run.py --doctor
python3 run.py rate --input counts.csv --field M_base --nt-mode per_tube
python3 run.py optimize --strain card.yaml      # protocol card optional
```

Windows: use `run.cmd` instead of `python3 run.py`.

## What it does

- `rate`: raw `r,Nt` counts → μ ± 95% CI (MSS-MLE, per-tube Nt, partial plating, Excel input)
- `compare`: two-batch rate comparison (LRT + fold-change CI)
- `power`: sample size / power (1σ or 95% CI conventions)
- `optimize`: three-knob window optimization with two hard protection caps (`P_mutator_max`, `P_protected_max`)
- `design-induction` / `fit-induction` / `fit-fc` / `fit-tauoff`: dose axis, dose-response and decay fitting

## Docs

Start with [`fluctuation-analysis/README.md`](fluctuation-analysis/README.md), then:

- `docs/tutorial.md`, `docs/glossary.md` — getting started
- `docs/knowledge-base.md` — principles and "seven whys"
- `docs/experiment-design.md` — design checklist and sample size
- `docs/mutator-window.md` — window model, protection caps, passenger budget
- `docs/viable-count-problems.md` — what to do when viable counts are unreliable
- `references/` — command manual, intake questions, guardrails, marker extrapolation

## Limitations

- Single-target Luria–Delbrück model; rRNA multi-copy / enzymatic / efflux resistance cannot be extrapolated by target size.
- Plasmid dynamics are simplified (`plasmid_loss_per_gen`), not a full copy-number simulation.
- Default protection caps (10% / 5%) are heuristics, not standards.
- Chinese documentation; Windows/macOS not tested on real hardware in this repo.

## Development

```bash
cd fluctuation-analysis
python3 selftest.py        # 21 checks
```
