+++
menuTitle = "Dlab Filesystem"
title = "Accessing the Dlab filesystem"
weight = 6
+++


Some Gradle tasks make use of the Dlab filesystem
to store files, like build artifacts.


## Configuring the Gradle build script

To configure Gradle to access these files, configure Gradle
by creating a `gradle.properties` file in the project root
with this content:

```properties
dlab.user = your-duke-cs-username
```

{{% notice note %}}
`gradle.properties` files always contain information that's
specific to the local computer or specific developer, so these files
are never committed to the Git repository.
{{% /notice %}}

The build script uses an SSH client to access the remote filesystem,
and the only supported authentication mechanism right now is public keys.

By default, the build script will look for your private key
at `$HOME/.ssh/id_rsa` and your public key at `$HOME/.ssh/id_rsa.pub`.
If your key is at a different location, add the paths to your `gradle.properties`:
```properties
dlab.key.private = /path/to/key
dlab.key.public = /path/to/key.pub
```

The SSH protocol has some ability to negotiate key types between
the client and server. For most users, this process will work correctly.
But sometimes key type negotiation can fail. If that happens, try
explicitly setting your key type in your `gradle.properties`:
```properties
dlab.key.type = key-type
```
For many keys (particularly older keys), the correct type is `ssh-rsa`.

The build script will also use your `known_hosts` file to authenticate
SSH servers. If your `known_hosts` file is not at `$HOME/.ssh/known_hosts`,
then add the path to your `gradle.properties`:
```properties
dlab.knownHosts = /path/to/known_hosts
```

## Usually it Just Works (in Linux)

For most users (particularly those using Linux), only the `dlab.user`
setting is required, and everything else will probably Just Work if your
SSH environment is configured in the usual way.

But if the build script can't establish the SSH connection, the
console output should hopefully give you clues where the errors
might be, and suggestions for how to fix them.


## OSPREY locations in the DLab filesystem

### The Release Archive {#release-archive}

**Location:** `/www/donaldlab/software/osprey/releases`

This folder houses all of the OSPREY release files. It is used by OSPREY build tools to store build artifacts.
It is also publicly available over HTTP(s) and is used to serve downloads for OSPREY releases to users.

**URL:**: https://www.cs.duke.edu/donaldlab/software/osprey/releases/

Take care when placing files here, since they will be visible to everyone on the internet.


### The documentatinon website {#docsite}

**Location:** `/www/donaldlab/software/osprey/docs`

This folder houses the live version of the documentation website and is publicly available over HTTP(s).

**URL:**: https://www.cs.duke.edu/donaldlab/software/osprey/docs/

Take care when placing files here, since they will be visible to everyone on the internet.
