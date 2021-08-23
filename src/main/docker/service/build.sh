#!/bin/sh

# this script requires sudo due to the way docker works!

# cd to this script's directory
cd `dirname "$0"`

# use the project root as the build context for docker
context=../../../..

# TODO: expose options for tag?
tag="osprey/service:latest"

docker build -f Dockerfile $context --tag $tag

# make a `usdo` (user do) command so we can do things as the logged in user
user=`logname`
usdo="sudo -u $user"

# make the build dir
build=$context/build/docker
$usdo mkdir -p $build

# save the docker image somewhere we can find it,
# but don't save the file with root permissions
target=$build/osprey-service.tar
docker save $tag > $target
chown $user: $target
