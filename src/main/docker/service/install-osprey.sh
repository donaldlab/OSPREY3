#!/bin/sh

home=/opt/osprey

# find all of the installed versions
versions=`find $home/versions/* -type d -name "v*"`
for version in $versions; do

  # get just the vXX.YY part of the version from the folder path
  vname=`basename $version`

  # derive the port number from the version
  # eg vXX.YY -> 8000 + int(XXYY)
  portstr=`echo ${vname#v} | sed -e "s/\.//g"`
  port=`expr $portstr + 8000`

  echo "installing osprey $vname on port $port"

  # get version-specific directories
  instdir=$home/versions/$vname
  srvdir=/etc/services.d/osprey-service-$vname

  # write the start script for this version
  mkdir -p $srvdir
cat << EOF > $srvdir/run
#!/bin/sh

cd /opt/osprey

# add java back to the path, since we somehow don't get it under s6
# thankfully, these path manipulations also seem to work when we execute with the www-data user
PATH=/opt/java/openjdk/bin:\$PATH

# run osprey service as the www-data user
s6-setuidgid www-data versions/$vname/bin/osprey-service.sh --port $port

EOF
  chmod +x $srvdir/run

  # copy the s6 finish script
  # NOTE: don't use symlinks here or s6 will be unhappy
  cp $home/finish.sh $srvdir/finish
  chmod +x $srvdir/finish

  # make sure the start script is runnable
  chmod +x $instdir/bin/osprey-service.sh

done
