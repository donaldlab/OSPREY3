
# Building the Osprey Service

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

Therefore, the build process for the service has different steps
for each layer of the service.


## Layer 1: Building the Java server

Run the `serviceTar` task in Gradle:
```shell
./gradlew serviceTar
```
This will build the service and place the distribution archive at
```
build/releases/osprey-service-$VERSON-tbz2
```
where `$VERSION` is the value of the `versionService` variable
in `build.gradle.kts`.


## Layer 2: Building the Docker multiplexer

### Download releases

First, gather all the versions of the service you want to bundle
into the Docker container by running the `downloadServiceReleases` Gradle task:
```shell
./gradlew downloadServiceReleases
```

This will download all built releases of the service into your
`build/releases` folder.

NOTE: If `versionService` in `build.gradle.kts` matches any of the previously-built
releases, that previously-built release will overwrite the releases of that
version you have built locally. If this is not the desired outcome,
bump your `versionService` to a new version number, or rebuild your
service after downloading relases.


### Build the Docker container

To configure which versions of the service are built into the Docker
container, edit `/src/main/docker/service/Dockerfile`. Add a build
command for each version to include like this:
```docker
ADD build/releases/osprey-service-$VERSION.tbz2 versions/v$VERSION
```
replacing `$VERSION` with the version you wish to include.
Then update the Docker container label to list all the included versions
by appending to the list of comma-delimited versions:
```docker
LABEL "osprey.service.versions"="[$OLD_VERSION[0],$OLD_VERSION[1],$NEW_VERSION]"
```
Here, `$OLD_VERSION[*]` are the old versions already in the file,
and `$NEW_VERSION` is the new version you wish to add.

For example, to include version `17.4` in the next build, when version
`16.2` is already present, find these lines:
```docker
ADD build/releases/osprey-service-16.2.tbz2 versions/v16.2
LABEL "osprey.service.versions"="[16.2]"
```
and change them to:
```docker
ADD build/releases/osprey-service-16.2.tbz2 versions/v16.2
ADD build/releases/osprey-service-17.4.tbz2 versions/v17.4
LABEL "osprey.service.versions"="[16.2,17.4]"
```

Then run the Docker container build script in Linux with `sudo`.
That's just how Docker works.
```shell
sudo src/main/docker/service/build.sh
```
When the script completes, the Docker container should appear at
`build/docker/osprey-service-docker-$VERSION.tar.bz2`
where `$VERSION` is the value of the `versionService` variable in
`build.gradle.kts`.


### Build the Docker release archive

Once the Docker container itself is ready, build the release archive
using the `serviceDockerTar` task in Gradle:
```shell
./gradlew serviceDockerTar
```
This will build the docker container and place the distribution archive at
```
build/releases/osprey-service-docker-$VERSION.tar
```
where `$VERSION` is the value of the `versionService` variable
in `build.gradle.kts`.


### TODO: publish docker release to dlab archive?
