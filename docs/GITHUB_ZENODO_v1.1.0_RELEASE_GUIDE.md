# GitHub → Zenodo release steps for v1.1.0

1. Replace the repository contents with the files in this package, or upload only the changed/new files.
2. Confirm that the default branch contains:
   - `CITATION.cff`
   - `.zenodo.json`
   - `LICENSE`
   - updated `README.md`
   - updated `CHANGELOG.md`
   - validation scripts/data/figures
3. Commit the changes to `main`.
4. On GitHub, open **Releases** → **Draft a new release**.
5. Create tag `v1.1.0` targeting `main`.
6. Release title: `Version 1.1.0 — Independent validation release`.
7. Paste the contents of `RELEASE_NOTES.md`.
8. Publish the GitHub release.
9. In Zenodo, verify that the GitHub integration is enabled for this repository. The release should be imported automatically.
10. Open the new Zenodo version, verify creators, MIT license, description, keywords, and files.
11. Record the new version-specific DOI. Continue using concept DOI `10.5281/zenodo.21128166` when referring to all software versions.
12. Update the manuscript Code Availability statement with the final version-specific DOI if the journal requires an exact archived version.
