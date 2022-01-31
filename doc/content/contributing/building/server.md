+++
menuTitle = "Building Osprey Server"
weight = 1
+++


# Building the Osprey Server

Run gradle task `serverRelease`:
```shell
./gradlew serverRelease
```

The built archive should show up in the folder
`build/releases/osprey-server-$OS-$VERSION.tbz2`,
where `$OS` desrcibes your operating system,
and `$VERSION` is the current version of Osprey.
