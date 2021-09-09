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

  # stop the daemon
  systemctl stop osprey
  systemctl disable osprey

  # uninstall the daemon
  rm "$systemddir/osprey.service"
  systemctl daemon-reload

elif [ -d "$launchddir" ]; then

  # nope nope nope, done messing around with OSX
  echo launchd detected, OSX platform not supported
  exit 1

  echo launchd detected, uninstalling for OSX

  # stop the daemon
  # launchd can automatically re-start stopped daemons, so disable before stopping
  launchctl disable osprey
  launchctl stop osprey

  # uninstall the daemon
  plistfile="$launchddir/osprey.plist"
  launchctl unload "$plistfile"
  rm "$plistfile"

fi

echo Removing docker image ...
docker image rm osprey/service:$version

echo Osprey service uninstalled!
