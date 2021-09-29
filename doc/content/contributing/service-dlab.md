
# Osprey Service Deployment at DLab

The public Osprey service is running on olympias.cs.duke.edu,
in a Docker container, via the `docker-machine` tools, since
Docker isn't natively supported on OSX.

https://github.com/docker/machine

https://docs.docker.com/machine/

Why not just use Linux to host our services? Good question.

We currently haven't provisioned any web-facing Linux machines
at the moment, but we do have one OSX one. So that's what we have for now.

Tragically, Docker doesn't really support OSX very well at all.
To actually run Docker on OSX, they made the `docker-machine` tool
that automatically provisions a Linux virtual machine with Docker
inside of it, since you really need Linux to run Docker.

And since OSX seems to be actively hostile to running system-level
daemons, all of the Dockers processes are currently running under
the user `jmartin`. If you need to take over maintenance of the
Osprey services on Olympias for some reason, you might be able to
use the `docker-machine` instance already installed under the `jmartin`
user, or you might need to install a new instance under your own user.
Or maybe you can figure out how to get OSX to run a proper system daemon
without requiring a user to be logged in
(`launchd` seems to have a ton of roadblocks there).
I don't really know.


## Managing Docker itself

### Checking the status of Docker

```shell
docker-machine ls
```

You should get a result like
```
NAME             ACTIVE   DRIVER       STATE     URL                         SWARM   DOCKER      ERRORS
osprey-service   -        virtualbox   Running   tcp://192.168.99.101:2376           v19.03.12
```

If not, you'll need to create a new VM.


### Creating a new docker machine for the Osprey Service

This part has already been done on Olympias,
but should you ever need to re-create the VM, run this:
```shell
docker-machine create --driver virtualbox osprey-service
```


### Starting/Stopping Docker itself

```shell
docker-machine start osprey-service
docmer-machine stop osprey-service
```


## Managing the Osprey service

### Getting access to Docker inside the VM

Before the `docker` command will even work, you'll have to
call `docker-machine` to install the connection parameters
into your shell's environment:

```shell
eval "$(docker-machine env osprey-service)"
```

### Check service status

```shell
docker ps
```

If the service is running, you should get a result like this:
```
CONTAINER ID   IMAGE                COMMAND   CREATED              STATUS              PORTS                      NAMES
abe880347ba6   osprey/service:0.3   "/init"   About a minute ago   Up About a minute   0.0.0.0:44342->44342/tcp   osprey-service
```

If not, you'll need to start the service, or even install it.

If the Docker container is running, you can test the public
availability of the service with a few different methods:
 * Run `curl --insecure https://localhost:44342/about` in your local shell.
 
   (The `--insecure` is needed because the service uses
   self-signed certs, which most TLS clients will complain about.)
 * Point your browser to https://olympias.cs.duke.edu:44342/about
 
   (and ask your browser to please not care about the self-signed certs.)  
 * Try to use the service from the Osprey GUI.


### Installing the Docker image

Get a copy of the Osprey service Docker release and untar it
to your favorite folder. Osprey releases can be found in the dlab archive
at `/usr/project/dlab/www/donaldlab/software/osprey/releases`.

Let's say you un-tar'd the release to the folder `$OSPREY`.

Also make note of the version number you're installing.
Let's say it's `$VERSION`.

eg, if you un-tar'd the file `osprey-service-docker-0.3.tar`,
then the version would be `0.3`.

Then install the docker image:
```shell
docker load < "$OSPREY/osprey-service-docker-$VERSION.tar.bz2"
```

If it worked, docker should show the image when running this command:
```shell
docker image ls
```
You should see something like:
```
REPOSITORY       TAG       IMAGE ID       CREATED       SIZE
osprey/service   0.3       077b07c945de   5 hours ago   734MB
```


### Starting Osprey service

Since we're running the osprey service in the maddness that is the
OSX environment, we'll need to do some extra work beyond what the
start script does to get things working. The start script was only
designed for proper Linux environments.

To start the service, first run the start script in the `$OSPREY` folder:
```shell
./osprey-service --version $VERSION > osprey-service.log &
```
where `$VERSION` matches the version given by `docker image ls`.

NOTE: Make sure the Docker daemon is running and accessible. If not, see above.

NOTE: On `olympias`, the `$OSPREY` folder is currently:
`/Users/jmartin/osprey/service-docker`

Use `&` to detatch the process from your shell and use `>` to save
a log of messages from the service.

TODO: redirect stderr AND stdout??

Once the service is started, configure the virtual machine to forward
ports so the service is visible to the internet:
```shell
VBoxManage controlvm "osprey-service" natpf1 "tcp-port44342,tcp,,44342,,44342";
```

The port forwarding must be done every time the service is started.


### Stopping Osprey service

To stop the service, run this command:
```shell
docker stop osprey-service
```
Easy peasy.


## Troubleshooting

### The Osprey Service in the Docker VM becomes unresponsive a few minutes after starting.
Executing `docker-machine ls` shows:
```
NAME             ACTIVE   DRIVER       STATE     URL                         SWARM   DOCKER    ERRORS
osprey-service   -        virtualbox   Running   tcp://192.168.99.101:2376           Unknown   Unable to query docker version: Get "https://192.168.99.101:2376/v1.15/version": remote error: tls: bad certificate
```

#### Possible Solution
The clock is not synchronized between the guest and host OSes.

Synchronize the clocks with:
```shell
docker-machine ssh osprey-service "sudo date -u $(date -u +%m%d%H%M%y)"
docker-machine regenerate-certs osprey-service
```

See https://github.com/docker/machine/issues/3845#issuecomment-516792766
for more info.


### docker-machine becomes unresponsive
Commands like `docker-machine ls` and `docker-machine stop $VM` hang forever.

#### Possible Solution
The VirtualBox VM that docker-machine manages has become unresponsive.
Try stopping it directly with VirtualBox tools:
```shell
VBoxManage list vms
VBoxManage list runningvms
VBoxManage controlvm osprey-service poweroff
```
See http://manpages.org/vboxmanage for more info.

Then maybe docker-machine will work again.


### docker-machine somehow becomes paused a few minutes after starting
`docker-machine ls` shows:
```
NAME             ACTIVE   DRIVER       STATE    URL   SWARM   DOCKER    ERRORS
osprey-service   -        virtualbox   Paused                 Unknown
```

Trying to resume the VM manually shows:
```shell
$ VBoxManage controlvm osprey-service resume
VBoxManage: error: VM is paused due to host power management
VBoxManage: error: Details: code VBOX_E_INVALID_VM_STATE (0x80bb0002), component ConsoleWrap, interface IConsole, callee nsISupports
VBoxManage: error: Context: "Resume()" at line 410 of file VBoxManageControlVM.cpp
```

This seems to be a very old bug in Virtualbox:
https://www.virtualbox.org/ticket/15378


#### Possible Solution
Looks like the VirtualBox VM is not able to recover after a guest OS
hibernation. The only known workaround is to stop the VM:
```
VBoxManage controlvm osprey-service poweroff
```
and then restart it.

Maybe the best solution is to auto-restart the VM after waking up
from hibernation, but I have no idea how to do that.

Perhaps the simpler solution is to just disable hibernation for now:
```shell
sudo pmset sleep 0
```
