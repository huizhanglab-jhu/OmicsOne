# Agent Skills

The `skills/` folder contains AI-agent instructions. These are not end-user
manuals and should stay short, operational, and specific to the repository.

## Included Skills

- `omicsone`: repository orientation and runtime selection.
- `omicsone-differential-analysis`: volcano, differential testing, gene mapping,
  and enrichment workflows.
- `omicsone-cnv-correlation`: CNV/RNA/protein Spearman correlation and CNV
  correlation figures.
- `omicsone-mutation-figures`: mutation heatmaps and mutation-type plots.
- `omicsone-boxplots`: normal-vs-tumor gene boxplots.
- `omicsone-pathway-scatter`: protein/phosphosite pathway scatter plots.

## Maintenance Rules

- Update skills when agent behavior should change.
- Update docs when human-facing usage, APIs, or outputs change.
- Keep trigger text in each skill's YAML `description`.
- Keep detailed tables and schemas in `docs/`, not in `SKILL.md`.
- Keep fragile workflow details in `references/` files when a skill becomes too
  long.

## Review Checklist

- Does the skill point to the correct source files?
- Does it tell the agent which runtime surface to use?
- Does it mention common pitfalls?
- Does it avoid duplicating long human documentation?
- Does it avoid stale hard-coded paths unless those paths are actual defaults?
