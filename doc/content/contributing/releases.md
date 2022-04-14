+++
menuTitle = "Releases"
title = "Releasing OSPREY"
weight = 5
+++


Releasing a new version of OSPREY involves several different steps.


## 1. Write the code

This is probably the most self-explanatory part, so there isn't much to say here.


## 2. Test the code

Hopfully the software we're about to release works correctly. There are many different philosophies for testing
software, but one of the most reliable methods is some kind of automated testing.

OSPREY's automated tests can be found in the `/src/test` folder.

Java code tends to use the popular [JUnit](https://junit.org/junit4/) framework (mostly version 4).

Kotlin code tends to use the [Kotest](https://kotest.io/) framework.

Your code IDE can help you run tests for specific parts of the codebase. Just open one of the testing source files
in your IDE and find the button that launches the selected tests. For example, in IntelliJ IDEA, you should see
little green arrows just to the left of the lines of test code.

Most IDEs have built-in support for JUnit. But for Kotest, you may need to install [an additional plugin][kotest-intellij].

[kotest-intellij]: https://kotest.io/docs/intellij/intellij-plugin.html


## 3. Increment the version number

Once you're happy with the state of the code and its status as a correctly working program,
bump the version number for OSPREY. The main version number for OSPREY is recorded in `/build.gradle.kts`
in a line that looks something like:
```kotlin
version = "4.0"
```

### Tips for picking version numbers

Version numbers in software usually take the form `$MAJOR.$MINOR.$REVISION` where `$MAJOR`, `$MINOR`, and `$REVISION`
are all small integers starting with `0`. When `$REVISION` is `0`, it is often omitted entirely.
Here are some example versions:

```
0.4
1.0
4.10.5
```

There are a few different philosophies for picking version numbers for software. Historically, OSPREY
has used one of the oldest ones. If stated explicitly, that philosophy might look something like:

> Increment the major version after we re-write most of the software and it looks totally different from before. \
> Increment the minor version after we add new features, or modify old ones. \
> Increment the revision number after we fix a bug and want to release the fix quickly.

This philosophy works rather well for projects without formal support contracts, service level agreements, or otherwise
paying customers. It's particularly well-suited to the more relaxed academic environment.

Other philosophies exist too, like the popular [SemVer](https://semver.org/). However, the focus of SemVer seems to be
to increase trust with downstream consumers of libraries by using the version system to signal when updates will cause
breaking changes. This is a fantastic idea in theory, but in practice, it's often hard to predict when changes in
the upstream software will cause breakages in downstream software. Getting these predictions right requires enormous
amounts of discipline on the part of the library developers. This is particularly hard to achieve in the academic
environment, when often the developers of the software are not always professional software engineers.
Even professional software engineers get SemVer wrong much more often than you'd like, so it's important to set
realistic expectations for OSPREY.

For those reasons, semantic versioning is a nice-to-have goal, but shouldn't be considered an essential part of
the release process for OSPREY.

When in doubt, the safest way to bump the version number is to just increment the minor version by one.
And remember that the next number after `9` is `10`.


### Commit your changes

Once you've found a new version number that brings you joy, commit your change to the git `main` branch,
or merge your commits into the `main` branch. Then [add a tag][git-tag] to your commit with the new joyful
version number. For example, your tag might be named `v4.1`.

[git-tag]: https://git-scm.com/book/en/v2/Git-Basics-Tagging


## 4. Push the commits to GitHub

Once your code is saved for all eternity in Git, the next step is to copy the code to the public
[Git repository for OSPREY on GitHub](https://github.com/donaldlab/OSPREY3) using the usual [`git push`][git-push]
mechanism.

[git-push]: https://git-scm.com/docs/git-push

The code will need to be publicly available for the next steps, which use [Continuous Integration (CI)][ci] services hosted
in [The Cloud][the-cloud]. If truth-be-told, OSPREY development doesn't actually follow the philosophies of
continuous integration very much at all. However, CI services are usually pretty good at doing cross-platform
builds, which OSPREY relies on to make life easier for developers that may only have access to the one platform
that happens to be running on their laptop.

[ci]: https://docs.microsoft.com/en-us/devops/develop/what-is-continuous-integration
[the-cloud]: https://azure.microsoft.com/en-us/overview/what-is-the-cloud/


## 5. Build the release using Azure Pipelines

Once your code is saved for all eternity in Git, the next step is to build OSPREY and create the release files.


### Log into Azure Pipelines

OSPREY uses [Azure Pipelines][azure-pipelines] to do cross-platform builds. To use it, you'll need to log in.

[azure-pipelines]: https://azure.microsoft.com/en-us/services/devops/pipelines/

First, log into the [DLab site on Azure](https://dev.azure.com/donaldlab). When prompted to enter your username,
enter your Duke NetID followed by `@duke.edu`. For example, `id123@duke.edu`.

{{% notice warning %}}
If you already have a Microsoft account, don't try to use it here. It won't work. You *must* use a Duke account.
{{% /notice %}}

You may need to request access from the Azure organization owner first though. If approval is needed,
the Azure website should walk you through those steps during the login process. If you still need assistance,
try asking the [Azure organization owner](https://dev.azure.com/donaldlab/_settings/organizationOverview) for help.

Once logged in, you should see the OSPREY project and its pipelines. Specifically, we'll be using the
[donaldlab.OSPREY4.release][azure-release-pipeline] pipeline.

[azure-release-pipeline]: https://dev.azure.com/donaldlab/osprey/_build?definitionId=6


### Install the Azure CLI

Before you can start the builds for OSPREY, you'll need to [install the Azure CLI][azure-cli-install].

[azure-cli-install]: https://docs.microsoft.com/en-us/cli/azure/install-azure-cli

Once installed, you can verify that the installation worked by running the `azureCheckCli` Gradle task.
If successful, you should see a message like:
```
Logged into Azure as: $YOUR_NETID@duke.edu
```
where `$YOUR_NETID` is your Duke NetID.

If you logged into the website in the previous steps, you should already be logged into the CLI too.
If for some reason the CLI thinks you're not logged in, then simply log in on the CLI directly with the shell command:
```
az login --allow-no-subscriptions
```

{{% notice tip %}}
The `--allow-no-subscriptions` flag is needed to prevent a warning when logging in. Something something SSO.
It's magic. Who really knows how it works.
{{% /notice %}}


### Actually start the build

Once the Azure CLI is installed and working, actually start the OSPREY builds with the `azureBuild` task.
This task uses the Azure CLI to start the build using your logged-in user.
Or, you could login to the Azure Devops website, click on OSPREY project, and manually run the [donaldlab.OSPREY4.release][azure-release-pipeline] pipeline.
These two options do the same thing.

After starting, you can watch the progress of the builds on the [donaldlab.OSPREY4.release][azure-release-pipeline] website.
The logs there are particularly helpful when troubleshooting problems.

The builds take 10-20 minutes to finish, so keep working dilligently on other tasks in parallel so as not to waste
precious time advancing towards your PhD. Or, you know, just [amuse yourself][reddit] until it's done.
I'm not the boss of you.

[reddit]: https://www.reddit.com/

When the build finishes, Azure well helpfully send you an email notification to your Duke email address.
It might be `$YOUR_NETID@duke.edu`, or `$YOUR.$NAME@duke.edu`. Who could know for sure. Once you get the email,
feel free to check the status of the build.

If you get some kind of red error indicator, then (at least part of) the build failed. Dive into the logs to
find out why. The reason might be something wrong with your code. Or the reason might be something wrong with the
pipeline itself. If you need to modify the pipeline itself for some reason, you can do that by editing the
`/release-pipelines.yml` file. That file is helpfully self-documenting, so it should help you get started
making changes to the pipeline itself if that's what you need.

On the other hand, if you got all green indicators, then the build succeeded and you can continue with the next steps!


## 6. Download the build artifacts

Once the build finishes successfully, run the `azureDownloadArtifacts` Gradle task. This will download all the
release files from the most recent build in Azure (hopefully the one you just ran) into your `/build/releases` folder.
After running this task, your local git clone will look like you built all the relases locally, even the releases
for different platforms.


## 7. Build the documentation

Run the `buildDocsRelease` Gradle task to build the documentation that goes with your new OSPREY release.


## 8. Archive the releases

Run the `archiveReleases` task to upload all the shiny new releases to the [DLab release archive][dlab-release-archive].
This will archive all the releases built by Azure, as well as the new documentation you just built.
This step will make all the new release files available for download by the public too, but they won't have links on
the website yet.

[dlab-release-archive]: {{< ref "/contributing/dlab-filesystem#release-archive" >}}

Which brings us to ...


## 9. Update the website

The last step is to update the OSPREY website. This is as simple as running a few more Gradle tasks (in order):
 * `downloadDocsReleases`
 * `buildWebsite`
 * `deployWebsite`

Be sure to read [updating the website][docdoc] before updating the website,
to get an overall understanding of how the OSPREY documentation system and the website tools work.

[docdoc]: {{< ref "/contributing/documentation#updating" >}}
