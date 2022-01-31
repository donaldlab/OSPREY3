#!/bin/sh

# this script requires sudo due to the way docker works!

# cd to this script's directory
cd `dirname "$0"`

# use the project root as the build context for docker
context=../../../../..

# pick a compressor for bzip2
if command -v lbzip2; then
  zip=lbzip2
  echo "lbzip2 compressor found"
elif command -v bzip2; then
  zip=bzip2
  echo "bzip2 compressor found"
  echo "It's fine"
  echo "But if you want to bzip2 faster, consider installing lbzip2"
else
  echo "No compressor for bzip2 found. Try installing lbzip2 or bzip2"
  exit 1
fi

# get the current osprey service version by parsing the gradle build code
# NOTE: don't try to run Gradle here... we shouldn't trust Gradle with sudo
version=` \
  cat "$context/buildSrc/src/main/kotlin/osprey/build/service.kt" \
  | grep -E "^\s*const val version = \"[^\"]+\"" \
  | sed -E "s/^\s*const val version = \"([^\"]+)\".*$/\\1/" \
`
if [ -z "$version" ]; then
  echo "Osprey service version not found, this is probably a bug in this build script"
  exit 1
fi

echo Building osprey service v$version Docker image ...

tag="osprey/service:$version"

docker build -f Dockerfile $context --tag $tag || exit 1

# make a `usdo` (user do) command so we can do things as the logged in user
user=`logname`
usdo="sudo -u $user"

# make the build dir
build=$context/build/docker
$usdo mkdir -p $build

# save the docker image somewhere we can find it,
# but don't save the file with root permissions
echo Exporting Docker image ...
target=$build/osprey-service-docker-$version.tar.bz2
docker save $tag | $zip > $target || exit 1
chown $user: $target
echo Done!
