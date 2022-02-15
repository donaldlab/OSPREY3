#!/bin/sh

# this script requires sudo due to the way docker works!

# cd to this script's directory
cd `dirname "$0"`

# use the project root as the build context for docker
context=../../../../..

# get the build folder, it should already exist
build=$context/build/docker

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

# get the current osprey service version by parsing the versions.txt file written by gradle
versionsfile="$build/versions.txt"
if [ ! -f "$versionsfile" ]; then
  echo "Versions file not found. Try running gradle prep task first."
  exit 1
fi
version=`sed -n '1p' < "$versionsfile"`
versions=`sed -n '2p' < "$versionsfile"`

echo Building osprey service v$version Docker image ...

tag="osprey/service:$version"

docker build -f Dockerfile $context --tag $tag --build-arg VERSIONS="$versions" || exit 1

# make a `usdo` (user do) command so we can do things as the logged in user
user=`logname`
usdo="sudo -u $user"

# save the docker image somewhere we can find it,
# but don't save the file with root permissions
echo Exporting Docker image ...
target=$build/osprey-service-docker-$version.tar.bz2
docker save $tag | $zip > $target || exit 1
chown $user: $target
echo Done!
