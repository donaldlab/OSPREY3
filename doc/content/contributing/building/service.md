+++
menuTitle = "Service"
title = "Building OSPREY Service"
weight = 3
+++


A fully-deployed instance of the Osprey Service has two different layers.

1. The Java HTTP server(s) each providing a version of the service itself
2. The Docker container with an HTTPs reverse proxy
   that multiplexes different versions of the service

This two-layer approach has a few benefits:

 * Web security (like HTTPs) is provided by extra software inside the
   Docker container and the Java code itself need not worry much
   about these details.
 * Multiple versions of the service can be hosted simultaneously,
   and clients can choose the version of the service they support.
   This is especially useful when inevitably older service clients
   will still expect to work with older versions of the service.
 * All of these features are bundled into a single Docker container,
   which can be easily deployed on any Docker-compatible host.
 
{{% notice note %}}
Due to the way Docker works, building OSPREY Service is only supported in Linux.
Additionally, you will need `sudo` access on the machine where you will be building the Docker container.
{{% /notice %}}

Therefore, the build process for the service has different steps
for each layer of the service.


## Layer 1: Building the Java server

Run the `serviceRelease` task in Gradle.
This will build the service and place the distribution archive at
```
build/releases/osprey-service-$VERSON.tbz2
```
where `$VERSION` is the value of the `BuildService.version` constant
in `buildSrc/src/main/kotlin/osprey/build/service.kt`.

{{% notice tip %}}
To see how to run tasks in Gradle, go to [Running Gradle Tasks]({{< ref "gradle#running-tasks" >}})
{{% /notice %}}


## Layer 2: Building the Docker multiplexer

### Preparation

First prepare the Docker build process by running the Gradle task `serviceDockerPrep`.
This will download all previously-built releases of the service into your
`build/releases` folder.

{{% notice note %}}
If the current `BuildService.version` matches any of the previously-built
releases, that previously-built release will overwrite the releases of that
version you have built locally. If this is not the desired outcome,
bump your `BuildService.version` to a new version number, or rebuild your
service after downloading relases.
{{% /notice %}}


### Build the Docker container

Then run the Docker container build script in Linux. You'll need to run the command with `sudo`.
That's just how Docker works.
```shell
sudo buildSrc/src/main/docker/service/build.sh
```
When the script completes, the Docker container should appear at
`build/docker/osprey-service-docker-$VERSION.tar.bz2`
where `$VERSION` is the value of the `BuildService.version`.

{{% notice note %}}
If the docker daemon is not started yet, start it with the command:\
`sudo systemctl start docker`
{{% /notice %}}


### Build the Docker release archive

Once the Docker container itself is ready, build the release archive
using the `serviceDockerRelease` task in Gradle.
This will build the docker container and place the distribution archive at
```
build/releases/osprey-service-docker-$VERSION.tar
```
where `$VERSION` is the value of the `BuildService.version`.


## Customizations to AmberTools

To remove some limitations in the [AmberTools][ambertools] used by OSPREY Service,
and make it more friendly to being called as a library, we have made some minor
customizations to the source code:

[ambertools]: https://ambermd.org/AmberTools.php

The patched AmberTools binaries are already included in the OSPREY git repository,
so applying the customizations are not necessary to build OSPREY.
However, if a developer needed to change the customizations for some reason,
here are instructions for patching AmberTools:

1. [Download AmberTools 19](https://ambermd.org/GetAmber.php) in source code format.
   Versions newer than 19 have not been tested.\
   {{% notice note %}}
   Tragically, the Amber group doesn't make it easy to download AmberTools19 anymore.
   A copy of the install files for AmberTools19 has been backed up on the
   [DLab filesystem]({{< ref "/contributing/dlab-filesystem" >}}) at `/Code/AmberTools19.tar.bz2`,
   in case we need it.
   {{% /notice %}}
2. Upack the AmberTools19.tar.bz2 file
3. export AMBERHOME=/path/to/unpacked/folder
4. `./configure --skip-python gnu`
    yes to all patches
5. apply patches in `progs/ambertools/patches` folder
6. `make install` (the build can be slow, e.g., `-j 4` can speed things up a lot)
7. copyable binaries are at `bin/` and `bin/to_be_dispatched`

More information about the patching process can be found in `progs/ambertools/patches/readme.txt`

To compile debug builds for testing:

1. `make clean`
2. `make AMBERBUILDFLAGS='-O0 -g' $MAKE_TARGET`

where `$MAKE_TARGET` is the Make target you're trying to build.
