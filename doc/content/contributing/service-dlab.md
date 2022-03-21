+++
menuTitle = "Dlab Osprey Service Deployment"
weight = 4
+++


# Osprey Service Deployment at DLab

The public Osprey service is running on `grisman-36.cs.duke.edu`,
in a [Docker][docker] container managed by a [systemd][systemd] service.

[docker]: https://www.docker.com/
[systemd]: https://systemd.io/


## Testing service availability

You can test the public availability of the service with a few different methods:


### At the command line

SSH into `grisman-36` and run:
```shell
curl --insecure https://localhost:44342/about
```

(The `--insecure` option is needed because the service uses
self-signed certificates, which most TLS clients will complain about.)

If the service is running, you should get a response like:
```
Osprey service running! Supported service versions: v1.0
```

### From a browser

Point your browser to https://grisman-36.cs.duke.edu:44342/about

(and ask your browser to please not care about the self-signed certs.)


### In OSPREY Desktop

Just use the features of OSPSREY Desktop and look for any errors.


## Check service status in systemd

OSPREY Service is installed as a systemd service named `osprey`.
Running:
```shell
systemctl status osprey
```
should give a result like this, if the service is running:
```
● osprey.service - Osprey Service
     Loaded: loaded (/lib/systemd/system/osprey.service; enabled; vendor preset: enabled)
     Active: active (running) since Wed 2022-03-16 20:22:35 UTC; 4 days ago
   Main PID: 43353 (osprey-service)
      Tasks: 23 (limit: 154591)
     Memory: 18.9M
     CGroup: /system.slice/osprey.service
             ├─43353 /bin/sh /opt/osprey-service/osprey-service --version 1.0
             └─43355 docker run --name osprey-service --rm -e S6_READ_ONLY_ROOT=1 --read-only --tmpfs /var:rw,ex>
```


### Check service status in Docker

To see if the OSPREY Service is running in Docker, run:
```shell
sudo docker ps
```

{{% notice note %}}
Interacting with the Docker daemon requires `sudo` privileges.
{{% /notice %}}

If the service is running, you should get a result like this:
```
CONTAINER ID   IMAGE                COMMAND   CREATED      STATUS      PORTS                                           NAMES
b2924c25f4de   osprey/service:1.0   "/init"   4 days ago   Up 4 days   0.0.0.0:44342->44342/tcp, :::44342->44342/tcp   osprey-service
```

If the service is not running, check to see that the Docker image for the service
is installed with:
```shell
sudo docker image ls
```

If the image is installed, you should see a response like:
```
REPOSITORY       TAG       IMAGE ID       CREATED       SIZE
osprey/service   1.0       2cd8e553d4b7   2 weeks ago   740MB
```
If the image is not installed, you will need to (re)install the OSPREY Service.


### Installing the Docker image

Downloads and installation instructions can be found on the [Main Page]({{< ref "/#osprey-service" >}}).


### Starting and Stopping Osprey service

Since OSPREY Service is managed by systemd in a service called `osprey`,
starting and stopping it is as simple as:
```shell
sudo systemctl start osprey
```
and
```shell
sudo systemctl stop osprey
```

By default, systemd is configured to start the OSPREY service at system boot time.
To change that behavior, use the `enable` and `disable` commands:
```shell
sudo systemctl enable osprey
```
and
```shell
sudo systemctl disable osprey
```
