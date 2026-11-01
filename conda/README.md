# Conda package

Build and test from the repository root:

```bash
conda create -n eac-build --override-channels -c conda-forge conda-build anaconda-client
conda activate eac-build
conda build conda/ --override-channels -c conda-forge --no-anaconda-upload --output-folder conda-dist
```

The recipe uses the versioned release source archive and its SHA256. Tests check
the installed CLI, offline review page and filtering workflow. The noarch package
supports Linux and macOS.

Publish the tested artifact:

```bash
anaconda upload --user westraingroup /absolute/path/to/eukasmchecker-VERSION-BUILD.conda
```

Keep credentials outside the repository. Verify installation from the public
channel after publishing.
