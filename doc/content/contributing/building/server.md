+++
menuTitle = "Server"
title = "Building OSPREY Server"
weight = 1
+++


## Prerequisites

Building OSPREY Server requires some extra Python development tools that are
not typically installed by default. Specifically, you'll need these Python packages:

 * [`wheel`](https://pypi.org/project/wheel/)
 * [`setuptools`](https://pypi.org/project/setuptools/)

Tragically, there's no universal way to install packages in Python. Python has somewhat of
a fragmented ecosystem. But generally the [`pip`][pip] tool works well, athough there are [many][conda] [other][pipx]
[alternatives][setuppy]. If you're a local administrator, you can try the simplest command:
```shell
pip install wheel setuptools
```
which might need a preceeding `sudo` on Linux. This is how [pypi.org](https://pypi.org) recommends
to install packages by default. If you're not a local administrator, you can do a
"user" installation using `pip`, which is good for shared computing environments like research clusters:
```shell
pip install --user wheel setuptools
```
Or if the `pip` command isn't on your path or is for the older version 2 of Python, you might try `pip3` instead.
Or use the Python module invocation instead of calling `pip` directly from the command line:
```shell
python -m pip install wheel setuptools
```
With that last one, feel free to add the `--user` flag if you need it. Or switch to the `python3`
executable if the default `python` executable ends up being Python 2 for some reason.
On Windows, the Python executable might even be `py`!

Easy, right?

[pip]: https://pip.pypa.io/en/stable/
[conda]: https://docs.conda.io/en/latest/
[pipx]: https://pypa.github.io/pipx/
[setuppy]: https://docs.python.org/3/distutils/setupscript.html

For more information on the unfolding saga of how to install Python packages in modern times,
[PyPI has a detailed guide](https://packaging.python.org/en/latest/tutorials/installing-packages/).


## Building

Once you're all set up, run gradle task `serverRelease`:
```shell
./gradlew serverRelease
```

The built archive should show up in the folder
`build/releases/osprey-server-$OS-$VERSION.tbz2`,
where `$OS` desrcibes your operating system,
and `$VERSION` is the current version of Osprey.
