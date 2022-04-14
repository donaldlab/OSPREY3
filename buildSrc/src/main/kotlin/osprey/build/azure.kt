package osprey.build

import org.gradle.api.Project
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.getValue
import org.json.JSONArray
import org.json.JSONObject
import java.io.ByteArrayOutputStream
import java.net.URL

import osprey.*


fun Project.makeAzureTasks() {

	fun String.quote(): String =
		"\"$this\""

	fun String.toJsonObject(): JSONObject =
		JSONObject(this)

	fun String.toJsonArray(): JSONArray =
		JSONArray(this)


	fun cli(cmd: List<String>): String {

		val cmdline = listOf("az") + cmd

		// log the command before running it
		println(cmdline.joinToString(" ") {
			if (" " in it) {
				it.quote()
			} else {
				it
			}
		})

		// run the command, capture stdout
		val out = ByteArrayOutputStream()
		exec {
			commandLine = cmdline
			setStandardOutput(out)
		}

		// parse the resulting JSON
		return out.toString(Charsets.UTF_8)
	}

	fun cli(vararg cmd: String) =
		cli(cmd.toList())


	val azureCheckCli by tasks.creating {
		group = "release"
		description = "Checks that the Azure CLI is available and a user is logged in"
		doLast {

			try {
				// try to get the logged in user
				val json = cli("ad", "signed-in-user", "show").toJsonObject()
				val user = json.getString("userPrincipalName")
				println("Logged into Azure as: ${user}")
				if (!user.endsWith("@duke.edu")) { throw Error("") }
			} catch (t: Throwable) {
				throw Error("""
					|Azure CLi check failed. Make sure that:
					|	1. Azure CLI is installed on your machine:
					|		https://docs.microsoft.com/en-us/cli/azure/install-azure-cli
					|	2. You are logged into your Duke account:
					|		Run `az login --allow-no-subscriptions`
					|		Login with your Duke NetID, for example: `abc123@duke.edu`
				""".trimMargin())
			}
		}
	}


	val azureOrganization = "https://dev.azure.com/donaldlab"
	val azureProject = "osprey"
	val azurePipeline = "donaldlab.OSPREY4.release"
	// TEMP
	//val azureBranch = "main"
	val azureBranch = "docsite-azure-pipelines"


	fun cliPipelines(cmd: List<String>): String {

		// split cmd into required and optional args
		val i = cmd.indexOfFirst { it.startsWith("--") }
		val cmdRequired: List<String>
		val cmdOptional: List<String>
		if (i >= 0) {
			cmdRequired = cmd.subList(0, i)
			cmdOptional = cmd.subList(i, cmd.size)
		} else {
			cmdRequired = cmd
			cmdOptional = emptyList()
		}

		// build the command
		return cli(listOf("pipelines")
			+ cmdRequired
			+ listOf(
				"--organization", azureOrganization,
				"--project", azureProject
			)
			+ cmdOptional
		)
	}

	fun cliPipelines(vararg cmd: String) =
		cliPipelines(cmd.toList())


	// command reference:
	// https://docs.microsoft.com/en-us/cli/azure/?view=azure-cli-latest
	// https://docs.microsoft.com/en-us/cli/azure/pipelines?view=azure-cli-latest

	tasks.create("azureBuild") {
		group = "release"
		description = "Starts the Azure Pipeline to build all the Osprey release artifacts"
		dependsOn(azureCheckCli)
		doLast {
			cliPipelines("run",
				"--name", azurePipeline,
				"--branch", azureBranch
			)
			println("""
				|
				|
				|Azure Pipelines build submitted!
				|
				|The build usually takes about 15 minutes to finish.
				|When it's done, you'll get an email with the results.
				|
				|You can also track the progress of the build here:
				|$azureOrganization/$azureProject/_build?definitionId=6
				|
			""".trimMargin())
		}
	}


	tasks.create("azureDownloadArtifacts") {
		group = "release"
		description = "Download the latest artifacts from the Azure Pipelines"
		doLast {

			// lookup the pipeline id
			val pipelineId = cliPipelines("list").toJsonArray()
				.map { it as JSONObject }
				.find { it.getString("name") == azurePipeline }
				?.getInt("id")
				?: throw Error("can't find pipeline: $azurePipeline")

			// get the latest run
			val runs = cliPipelines("runs", "list",
				"--pipeline-ids", pipelineId.toString(),
				"--branch", azureBranch,
				"--query-order", "FinishTimeDesc",
				"--top", "1"
			).toJsonArray()
			val run = runs
				.takeIf { it.length() > 0 }
				?.getJSONObject(0)
				?: throw Error("no pipeline runs from which to download artifacts")

			val buildId = run.getInt("id")

			// make sure it's finished successfully
			if (run.getString("status") != "completed") {
				throw Error("pipeline run not completed yet")
			}
			if (run.getString("result") != "succeeded") {
				throw Error("""
					|pipeline run was not successful. Fix the problem and then try again:
					|$azureOrganization/$azureProject/_build/results?buildId=$buildId&view=results
				""".trimMargin())
			}

			// get the artifact groups
			val artifactGroups = cliPipelines("runs", "artifact", "list",
				"--run-id", buildId.toString()
			).toJsonArray()

			releasesDir.createFolderIfNeeded()

			for (artifactGroup in artifactGroups.filterIsInstance<JSONObject>()) {
				val name = artifactGroup.getString("name")

				// download the zip file of the artifacts in this group
				val zipPath = releasesDir.resolve("$name.zip")
				println("downloading $zipPath ...")
				URL(artifactGroup.getJSONObject("resource").getString("downloadUrl"))
					.download(zipPath)

				// unzip the file (and flatten the folders)
				println("unpacking ${zipPath.fileName} ...")
				copy {
					from(zipTree(zipPath).files)
					into(releasesDir)
				}
				zipPath.deleteFile()
			}

			// the linux desktop build comes out with a weird name (thanks, debian)
			// eg: osprey-desktop_4.0-1_amd64.deb
			// jpackage seems to hard-code this format with no configuration options:
			// https://github.com/openjdk/jdk/blob/master/src/jdk.jpackage/linux/classes/jdk/jpackage/internal/LinuxDebBundler.java#L118
			// the eg `-1` bit is the release number
			// since we're packaging our own releases, we don't really need it
			// so rename the file to eg: osprey-desktop-4.0-amd64.deb
			val linuxDebPattern = Regex("^([^_]+)_([0-9.]+)-[0-9]+_([^.]+)\\.deb$")
			for (file in releasesDir.listFiles()) {
				val match = linuxDebPattern.matchEntire(file.fileName.toString()) ?: continue
				val pack = match.groups.get(1)!!.value
				val version = match.groups.get(2)!!.value
				val arch = match.groups.get(3)!!.value
				file.rename("$pack-$version-$arch.deb")
			}
		}
	}
}
