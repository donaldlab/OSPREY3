#!/bin/sh

# this script installs the Osprey service docker image and starts it
# it must be run with sudo

# get the install dir (this script dir) as an absolute path
here=`realpath "$0" | xargs dirname`

# get the service version from the docker image filename
version=` \
  ls "$here/osprey-service-docker-"* \
  | xargs -I "{}" basename "{}" ".tar.bz2" \
  | sed -E "s/^osprey-service-docker-(\\.*)/\\1/" \
`

# import the image into docker
image="$here/osprey-service-docker-$version.tar.bz2"
echo Importing image $image into Docker ...
docker load < "$image" || exit 1

# install service to the OS's init system
systemddir="/usr/lib/systemd/system"
launchddir="/Library/LaunchDaemons"
if [ -d "$systemddir" ]; then

  echo systemd detected, installing for linux

  # make the unit file
  # https://www.freedesktop.org/software/systemd/man/systemd.unit.html
  cat << EOF > "$systemddir/osprey.service" || exit 1
[Unit]
Description=Osprey Service

[Service]
Type=simple
Restart=always
RestartSec=1
WorkingDirectory=$here
ExecStart=$here/osprey-service --version $version

[Install]
WantedBy = multi-user.target
EOF

  # poke systemd to recognize the new daemon and start it
  systemctl daemon-reload
  systemctl start osprey

  # start the daemon at boot too
  systemctl enable osprey


elif [ -d "$launchddir" ]; then

  echo launchd detected, installing for OSX

  # https://developer.apple.com/library/archive/documentation/MacOSX/Conceptual/BPSystemStartup/Chapters/CreatingLaunchdJobs.html

  # TODO

fi

echo Osprey service installed!
echo Try pointing a browser to https://localhost:44342/about to check
