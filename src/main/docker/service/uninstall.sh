#!/bin/sh

# get the install dir (this script dir) as an absolute path
here=`realpath "$0" | xargs dirname`

# get the service version from the docker image filename
version=` \
  ls "$here/osprey-service-docker-"* \
  | xargs -I "{}" basename "{}" ".tar.bz2" \
  | sed -E "s/^osprey-service-docker-(\\.*)/\\1/" \
`

# uninstall the service from the OS's init system
systemddir="/usr/lib/systemd/system"
launchddir="/Library/LaunchDaemons"
if [ -d "$systemddir" ]; then

  echo systemd detected, uninstalling for linux

  systemctl stop osprey
  systemctl disable osprey
  rm "$systemddir/osprey.service"
  systemctl daemon-reload

elif [ -d "$launchddir" ]; then

  echo launchd detected, uninstalling for OSX

  # TODO

fi

echo Removing docker image ...
docker image rm osprey/service:$version

echo Osprey service uninstalled!
