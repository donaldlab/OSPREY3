+++
title = "OSPREY Service on Linux"
weight = 8
+++


OSPREY Service is designed to run on Linux.
It is not supported on any other operating system.

OSPREY Service has been tested on Ubuntu Linux,
but it likely will work with any Linux distribution that uses
systemd, including distributions based on Debian or Fedora.

Since OSPREY Service uses Docker, it may even work on other Linux
distributions not based on Debian or Fedora, but you will have to
set up automatic boot-time startup using your distribution's init system.
The install script will fail if systemd is not available.


## 1. Prerequisites

OSPREY Service is designed to run with [Docker][docker]

[docker]: https://www.docker.com/

First, [install Docker for Ubuntu Linux][docker-ubuntu].

[docker-ubuntu]: https://docs.docker.com/engine/install/ubuntu/


## 2. Download OSPREY Service release

If you haven't already, download the OSPREY Service release.

[OSPREY downloads](/install/versions)

Extract the release `.tar` file you downloaded into your favorite folder.
Let's call that folder `$DIR`.

```shell
cd $DIR
tar -xvf osprey-service-docker-$VERSION.tar
```
where `$VERSION` is the version of OSPREY service that you downloaded.

Then, delete the `.tar` file we downloaded. It's not needed anymore.
```shell
rm osprey-service-docker-$VERSION.tar
```
It's actually very important to delete this file. Otherwise,
the installation script gets confused and will fail.


## 3. Install OSPREY Service

Run the install script in `$DIR` with `sudo`.
```shell
sudo ./install.sh
```

Importing the image into Docker can take a few minutes to finish.
Then the rest of the installation should finish very quickly.

When the script finishes, you can test the installation by running this command:
```shell
curl -k https://localhost:44342/about
```

You should see a response like:
```
Osprey service running! Supported service versions: $VERSIONS
```

## 4. Managing OSPREY Service

The install script created a systemd service called `osprey` to run
OSPREY automatically at boot time. You can manage this service using
the usual systemd commands.

To check on the status of the OSPREY service, run:
```shell
systemctl status osprey
```

If OSPREY is running, you should see a response like:
```
● osprey.service - Osprey Service
     Loaded: loaded (/lib/systemd/system/osprey.service; enabled; vendor preset: enabled)
     Active: active (running) since Wed 2022-03-16 19:55:35 UTC; 16min ago
   Main PID: 42245 (osprey-service)
      Tasks: 19 (limit: 154591)
     Memory: 18.8M
     CGroup: /system.slice/osprey.service
             ├─42245 /bin/sh /opt/osprey-service/osprey-service --version 1.0
             └─42247 docker run --name osprey-service --rm -e S6_READ_ONLY_ROOT=1 --read-only --tmpfs /var:rw,exec --mount type=bind,src=/tmp,dst=/tmp --vol>
```

To stop the OSPREY service, run:
```shell
sudo systemctl stop osprey
```

To start the OSPREY service again, run:
```shell
sudo systemctl start osprey
```

To prevent the OSPREY service from starting at boot, run:
```shell
sudo systemctl disable osprey
```

To start OSPREY at boot time again, run:
```shell
sudo systemctl enable osprey
```
