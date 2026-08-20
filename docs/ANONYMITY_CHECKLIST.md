# Anonymity checklist

This repository will be linked from a double-blind ICLR 2027 submission. Work
through every item before making it public. `make anonymity-check` automates the
mechanical greps; the judgement calls below are not automatable.

Run the check well before the deadline. Several items — git history in
particular — require rewriting history, which is disruptive if discovered late.

---

## Git

- [ ] **Author strings.** `git log --format='%an <%ae>' | sort -u` shows no real
      name, no institutional email, no username that identifies a person. This
      includes committer as well as author fields.
- [ ] **Commit messages** mention no person, lab, institution, cluster, grant, or
      internal ticket system.
- [ ] **Branch and tag names** carry no identifying strings.
- [ ] **Remotes.** No institutional git host in `.git/config` in the published
      copy.
- [ ] If any of the above fails, history is rewritten (`git filter-repo`) and the
      repository force-pushed to a clean location — not merely amended on the tip
      commit, which leaves the rest of the history intact.

Set the anonymous identity **before** the first commit:

```bash
git config user.name "Anonymous" && git config user.email "anonymous@example.com"
```

---

## Files

- [ ] **No author names** anywhere in tracked files — source, configs, docs,
      notebooks, figures.
- [ ] **No institution names**, including in acknowledgements, funding notes, or
      dataset access notes ("obtained via the X consortium").
- [ ] **No absolute paths** revealing a home directory, username, or cluster:
      `git grep -nE '/(home|Users|scratch|nfs|mnt)/[A-Za-z0-9._-]+'`
- [ ] **No cluster or machine names** in configs, scripts, SLURM files, or
      logging output.
- [ ] **LICENSE** copyright holder reads `Anonymous`.
- [ ] **`pyproject.toml`** `authors` field reads `Anonymous`; `Repository` URL
      points at an anonymised host.
- [ ] **No internal URLs** — wikis, issue trackers, shared drives, dashboards.

---

## Notebooks

- [ ] **All output cells cleared.** Outputs leak paths, hostnames, usernames in
      tracebacks, and figure metadata. The `nbstripout` pre-commit hook enforces
      this; verify it is actually installed (`pre-commit install`), since a hook
      that was never installed silently enforces nothing.
- [ ] **Execution counts cleared.**
- [ ] **No kernel metadata** naming a personal environment or conda prefix.

---

## Figures and results

- [ ] **No author or institution text** in any figure, including embedded font
      metadata and PDF producer strings.
- [ ] **Run manifests** (`results/*/manifest.json`) contain no username,
      hostname, home directory, or cluster name. `utils/io.py::write_manifest`
      is responsible for this; check the actual output rather than trusting the
      docstring.
- [ ] **No results committed** beyond `results/.gitkeep`.

---

## Prior work

This is the item most likely to be missed, because it looks like good scholarly
practice rather than a leak.

- [ ] **Related prior work by the same group is cited in third person**, exactly
      as any other prior work would be. In particular `b-LOAD` (CIKM 2026),
      whose Meek/MPDAG closure is this project's starting point.
- [ ] **No "our previous work", "we showed in", "our earlier paper"** anywhere in
      the repository or the paper.
- [ ] **No vendored code** from a group repository carrying its origin in
      comments, headers, or module docstrings. If the b-LOAD closure is used as
      a starting point, it is reimplemented and attributed as a citation, not
      copied with provenance comments.
- [ ] **`docs/`** does not describe the lineage of the work in a way that
      identifies the authors.

---

## Dependencies

- [ ] **No private package indexes** or institutional wheels in
      `pyproject.toml`, `requirements.txt`, or any lockfile.
- [ ] **No git dependencies** pointing at an institutional host.

---

## Final pass

- [ ] `make anonymity-check` is clean.
- [ ] A fresh clone into a new directory installs and passes tests without any
      path outside the repository.
- [ ] Someone other than the author reads the README as if they had never seen
      the project, and reports nothing identifying.
