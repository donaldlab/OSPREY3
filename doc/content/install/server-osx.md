+++
title = "OSPREY Server on OSX"
weight = 6
+++

## 1. Prerequisites

First, install [Python][python] and [pip][python-pip], if not installed already.

[python]: https://www.python.org/
[python-pip]: https://pip.pypa.io/en/stable/

Version 3 or higher must be used. Make sure to get the 64-bit verison,
and not a 32-bit version.

You can download Python from the official Python website, or you can
install Python via [Homebrew][homebrew]:
```shell
brew install python
```

[homebrew]: https://brew.sh/


## 2. Download OSPREY Server release

If you haven't already, download the OSPREY Server release.

[OSPREY downloads]({{< ref "/install/versions" >}})

Extract the release `.tar.bz2` file you downloaded into your favorite folder.
From that folder, run the installation script to install OSPREY:
```shell
./install.sh
```


## 3. Test your installation

Start an interactive Python session, then run:
```python
import osprey
osprey.start()
```

If successful, should should be greeted with a message something like the following:
```
OSPREY 3.2, Python 3.6.9, Java 11.0.6
Using up to 1024 MiB heap memory: 128 MiB for garbage, 896 MiB for storage
```


## 4. Updating to a new version of OSPREY Server

If you're updating from an older installation of Osprey,
just run the install script. There's no need
to explicitly uninstall the older version.


## 5. To Uninstall OSPREY Server

Just run the uninstall script:
```shell
./uninstall.sh
```
