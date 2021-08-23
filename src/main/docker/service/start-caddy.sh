#!/usr/bin/with-contenv sh
# NOTE: the above shebang is specific to s6
# it lets this service see all the environment variables set by the osprey-service launcher script (via docker)

cd /opt/osprey

# use s6 to run things as the user, like sudo, but user-do
usdo="s6-setuidgid www-data"

# configure folders for caddy
caddydir=/home/osprey/caddy
$usdo mkdir -p $caddydir/config
$usdo mkdir -p $caddydir/data
$usdo mkdir -p $caddydir/home
export XDG_CONFIG_HOME=$caddydir/config
export XDG_DATA_HOME=$caddydir/data
export HOME=$caddydir/home

# generate the SSL certificate for Caddy if needed
# this certificate will be self-signed, but that's ok, since we can make all the clients accept them anyway
# browsers will complain about self-signed certificates, but who cares... this isn't a browser app
keypath=$HOME/osprey-service
if [ ! -f "$keypath.key" ]; then
  $usdo openssl req -new -newkey rsa:4096 -x509 -sha256 -days 1 -nodes -out $keypath.crt -keyout $keypath.key << EOF || exit 1
.
.
.
Osprey
.
osprey.service
osprey@cs.duke.edu

EOF
  $usdo chmod go-rwx $keypath.key
  $usdo chmod go+r $keypath.crt
fi

# find all of the installed osprey versions
proxies=
ospreydir=/opt/osprey
versions=`find $ospreydir/versions/* -type d -name "v*"`
for version in $versions; do

  # get just the vXX.YY part of the version from the folder path
  vname=`basename $version`

  # derive the port number from the version
  # eg vXX.YY -> 8000 + int(XXYY)
  portstr=`echo ${vname#v} | sed -e "s/\.//g"`
  port=`expr $portstr + 8000`

  echo "proxying osprey $vname on port $port"

  # build caddy directives to set up a reverse proxy
  proxies="$proxies
  handle_path /$vname/* {
    reverse_proxy localhost:$port
  }"

done

# package the caddy directives into an environment variable
export CADDY_PROXIES="$proxies"

# NOTE: HTTPS_PORT should be set (via docker, then s6) in the osprey-service launcher
if [ $HTTPS_PORT -gt 0 ]; then
  echo "Starting Osprey service proxy on HTTPs port $HTTPS_PORT"
else
  echo "ERROR: No HTTPs port configured! $HTTPS_PORT"
  exit 1
fi

# run caddy as the www-data user
$usdo caddy run --config Caddyfile
