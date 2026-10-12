# Publishing RHZLS manually

The local repository has a `main` branch, committed source files, and a remote
named `origin` pointing to `https://github.com/MarcoSutti/RHZLS.git`.
The publication is performed by Marco Sutti.

## Create the remote repository

1. Sign in to GitHub as MarcoSutti and open https://github.com/new.
2. Choose repository name **RHZLS** and visibility **Public**.
3. Leave the initial README, licence and gitignore options unchecked: the local
   repository already provides these files.
4. Create the repository.

## Push the prepared history

From a PowerShell terminal in the RHZLS directory:

```powershell
git push -u origin main
git push origin v0.1.0
```

The first push may ask you to authenticate. Authenticate directly with GitHub;
do not place passwords or access tokens in source files.

If Git reports a directory-ownership warning in the original Codex workspace,
trust only this intended repository for the individual command:

```powershell
$repositoryPath = (Get-Location).Path
git -c "safe.directory=$repositoryPath" push -u origin main
git -c "safe.directory=$repositoryPath" push origin v0.1.0
```

## Publish the versioned release

Open the repository's **Releases** page, choose **Draft a new release**, select
the existing `v0.1.0` tag, and use the title **RHZLS 0.1.0**. Summarize the changes
from [CHANGELOG.md](../CHANGELOG.md). The prepared source ZIP can be attached as
an asset; Manopt is installed separately at the fixed tested version 8.0.

The source checkout includes nine reference PDFs. It excludes Manopt sources,
generated MAT files, local logs, private development archives and Git metadata.
Users reproduce all MAT results and figures with `run_all_experiments` after
selecting an external Manopt installation.

## Check automated tests and citation

The GitHub Actions workflow runs the regression suite on push and pull request.
A manual workflow input also runs all 58 experiments. Hosted CI results must be
checked after publication; local verification does not imply a hosted CI pass.

`CITATION.cff` identifies Marco Sutti as the software author and the Sutti and
Vandereycken 2021 paper as the preferred research citation. Check GitHub's
**Cite this repository** menu after the push.

For an archival DOI, optionally enable the repository in Zenodo before
publishing a release: https://help.zenodo.org/docs/github/. Add the assigned DOI
to the citation metadata after it exists; do not invent a DOI or replace an
already published version tag to change its history.
