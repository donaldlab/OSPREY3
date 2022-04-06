+++
menuTitle = "Desktop"
title = "Building OSPREY Desktop"
weight = 2
+++


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

{{% notice note %}}
The `.deb` packaging tools will only work in a Debian-based linux,
and not other distributions like Fedora/RHEL/CentOS/Rocky/Alma or Arch/Manjaro.
{{% /notice %}}

<!-- TODO: build RPMs for Fedora/RHEL/CentOS/Rocky/Alma/etc ? -->
<!-- TODO: find a way to build for Arch/Manjaro? ? -->


## for OSX:

OSPREY Desktop builds as a `.dmg` application installer in OSX.


## for Windows:

Before building, install the [WiX tools][wix] from Microsoft.

[wix]: https://wixtoolset.org/

Under Windows, OSPREY Desktop builds as an MSI installer.
