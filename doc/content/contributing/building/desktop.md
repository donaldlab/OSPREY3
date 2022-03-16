+++
menuTitle = "Building OSPREY Desktop"
weight = 2
+++


# Building OSPREY Desktop

Run gradle task `desktopRelease`:
```shell
./gradlew desktopRelease
```

The built archive should show up in the folder
`build/releases/osprey-desktop-$OS-$VERSION.tbz2`,
where `$OS` desrcibes your operating system,
and `$VERSION` is the current version of Osprey.


## for Linux:

OSPREY Desktop builds as a `.deb` package for Debian-based Linux
distributions, including Debian, Ubuntu, Mint, and others.

<!-- TODO: build RPMs for Fedora/RHEL/CentOS/Rocky/Alma/etc ? -->


## for OSX:

TODO: what does jpackage do on OSX?


## for Windows:

Before building, install the [WiX tools][wix] from Microsoft.

[wix]: https://wixtoolset.org/

Under Windows, OSPREY Desktop builds as an MSI installer.
