#!/bin/sh

# this script requires sudo due to the way docker works!

# cd to this script's directory
cd `dirname "$0"`

# use the project root as the build context for docker
context=../../../..

# pick a compressor for bzip2
if [ -x "`which lbzip2`" ]; then
  zip=lbzip2
  echo "lbzip2 compressor found"
elif [ -x "`which bzip2`" ]; then
  zip=bzip2
  echo "bzip2 compressor found"
  echo "It's fine"
  echo "But if you want to bzip faster, consider installing lbzip2"
else
  echo "No compressor for bzip2 found. Try installing lbzip2 or bzip2"
  exit 1
fi

# get the current osprey service version by parsing build.gradle.kts
# NOTE: don't try to run Gradle here... we shouldn't trust Gradle with sudo
version=` \
  cat "$context/build.gradle.kts" \
  | grep -E "^val versionService = \"[^\"]+\"" \
  | sed -E "s/^val versionService = \"([^\"]+)\".*$/\\1/" \
`
echo Building osprey service v$version Docker image ...

tag="osprey/service:$version"

docker build -f Dockerfile $context --tag $tag

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
docker save $tag | $zip > $target
chown $user: $target
echo Done!
