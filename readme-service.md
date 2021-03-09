
# Osprey Service

Provides functions for the cross-platform Osprey GUI that can't be run
directly on client computers, usually due to unsupported operating systems.


## Installation

### Dependencies

Osprey Service needs a Java runtime (v8+) and a fortran runtime.

To install on Ubuntu (or other Debian-based linux):
```shell
sudo apt install openjdk-11-jre-headless, libgfortran-7-dev
```

### Deploy application files

Choose a path to deploy the files. Let's refer to that path as `$DIR`.
Then extract the distribution archive file into `$DIR`.

Then make sure the `logs` directory is writeable by the `www-data` user:
```shell
cd $DIR
mkdir logs
sudo chown www-data:www-data logs
```


### Create a service on systemd

```shell
cd /etc/systemd/system
sudo vim osprey.service
```

Add this text to the file:
```
[Unit]
Description=Service for Osprey GUI
After=network.target

[Service]
Type=simple
Restart=always
RestartSec=1
User=www-data
WorkingDirectory=$DIR
ExecStart=$DIR/bin/osprey-service

[Install]
WantedBy = multi-user.target
```
Don't forget to replace `$DIR` in this file with the path you chose.

Then start the service:
```shell
sudo systemctl start osprey
```

Then configure the service to start a boot
```shell
sudo systemctl enable osprey
```

If you change the service file, get systemd to reload it:
```shell
sudo systemctl daemon-reload
```