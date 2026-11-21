# Using the EXISTING university repository

Do **not** run `git init` or make a new remote if your supervisor already provided a repository. Clone the existing remote and preserve its existing README, license, and branches.

## First-time setup (Windows PowerShell or Git Bash)

```bash
git clone <UNIVERSITY_REPOSITORY_URL>
cd <YOUR_CLONED_REPO_FOLDER>
git status
git branch -a
```

Review the repository. If the university requires work in a specific directory or default branch, use that layout. Otherwise, extract the **contents** of this starter ZIP into the cloned repository root. If a file already exists (especially README.md or .gitignore), manually merge its content instead of overwriting it.

```bash
git switch -c feature/step-01-project-scope
git status
git add docs/ config/ data/metadata/ notebooks/ src/ tests/ results/ thesis/ environment.yml
# Review staged files and add README/.gitignore only after merging existing versions
git diff --cached --stat
git diff --cached
git commit -m "docs: initialize landslide risk study and research scope"
git push -u origin feature/step-01-project-scope
```

Create a Pull Request / Merge Request into the repository's default branch if this is the agreed workflow. If you have no push rights, ask the repository owner which contribution process to use. If authentication is needed, follow the university's SSO/SSH/PAT procedures; never put a password or token in code.

## Subsequent steps
- Start from the accepted/default branch and create a short-lived branch such as `feature/step-03-study-area`.
- Make one coherent change, document inputs and outputs, run checks, inspect `git status` and `git diff --cached`.
- Commit to the branch and push; request review if needed.
- Avoid committing raw GIS rasters, commercial/restricted data, caches or secrets. `.gitignore` protects only files that are not already tracked.
- If data files really need to be shared, first confirm license/permissions and consider university storage or Git LFS only where approved and supported.
- Use repository tags or release notes to mark dissertation result snapshots, and capture package versions from a known working environment.

## Why the project repository is separate from CLIMADA
CLIMADA is a dependency; this repository holds **our analysis**, not copies of `climada_python` and `climada_petals`. Installation information is in `environment.yml` and README.

## Git basics after work is finished

```bash
git status
git add docs/ src/ notebooks/ config/ data/metadata/ results/figures/ results/tables/
git diff --cached --stat
git commit -m "feat: add study-area mapping workflow"
git push
```

Do not run `git add -A` blindly when working in directories with unreviewed data files.
