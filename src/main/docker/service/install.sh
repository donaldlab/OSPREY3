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
  servicefile="$systemddir/osprey.service"
  cat << EOF > "$servicefile" || exit 1
[Unit]
Description=Osprey Service
After=docker.service

[Service]
Type=simple
Restart=always
RestartSec=1
WorkingDirectory=$here
ExecStart=$here/osprey-service --version $version

[Install]
WantedBy = multi-user.target
EOF
  chmod go-w "$servicefile"

  # poke systemd to recognize the new daemon and start it
  systemctl daemon-reload
  systemctl start osprey

  # start the daemon at boot too
  systemctl enable osprey


elif [ -d "$launchddir" ]; then

  # this is too much trouble to get working now
  # maybe try again some other day ...
  echo launchd detected, OSX platform not supported
  exit 1

  echo launchd detected, installing for OSX
  # https://developer.apple.com/library/archive/documentation/MacOSX/Conceptual/BPSystemStartup/Chapters/CreatingLaunchdJobs.html

  # make the properties file
  plistfile="$launchddir/osprey.plist"
  cat << EOF > "$plistfile" || exit 1
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
  <dict>
    <key>Label</key>
    <string>osprey</string>
    <key>ProgramArguments</key>
    <array>
      <string>$here/osprey-service</string>
      <string>--version</string>
      <string>$version</string>
    </array>
    <key>KeepAlive</key>
    <true/>
    <key>UserName</key>
    <string>nobody</string>
  </dict>
</plist>
EOF
  chmod go-w "$plistfile"
  # tragically, launchd doesn't seem to support dependencies for always-on daemons
  # there doesn't seem to be any way to guarantee that docker is running first =(
  # does apple even WANT devs to make services for OSX servers?!?

  # install the daemon
  launchctl load "$plistfile"
  launchctl enable osprey

  # start the daemon
  launchctl start osprey

fi

echo Osprey service installed!
echo Try pointing a browser to https://localhost:44342/about to check
