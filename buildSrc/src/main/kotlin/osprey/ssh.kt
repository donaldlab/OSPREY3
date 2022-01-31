package osprey

import org.gradle.api.Project
import java.nio.file.Files
import java.nio.file.Paths
import com.jcraft.jsch.JSch
import com.jcraft.jsch.ChannelSftp
import com.jcraft.jsch.Logger as JschLogger
import com.jcraft.jsch.SftpProgressMonitor


/** configure and open an SFTP connection over SSH */
fun <T> Project.sftp(block: ChannelSftp.() -> T): T {

	// get user auth info
	val sshDir = Paths.get(System.getProperty("user.home")).resolve(".ssh")
	val user = project.propertyOrNull("dlab.user") as String?
		?: throw Error("no user configured. set `dlab.user = your-user` in gradle.properties")
	val keypriv = (project.propertyOrNull("dlab.key.private") as String?)
		?.let { Paths.get(it) }
		?: sshDir.resolve("id_rsa")
	val keypub = (project.propertyOrNull("dlab.key.public") as String?)
		?.let { Paths.get(it) }
		?: Paths.get("$keypriv.pub")
	val keytype = (project.propertyOrNull("dlab.key.type") as String?)
	val knownHosts = (project.propertyOrNull("dlab.knownHosts") as String?)
		?.let { Paths.get(it) }
		?: sshDir.resolve("known_hosts")

	// configure host info
	val host = "login.cs.duke.edu"
	val port = 22

	// API docs: http://www.jcraft.com/jsch/
	val jsch = JSch()
	jsch.addIdentity(keypriv.toString(), keypub.toString(), null)
	jsch.setKnownHosts(knownHosts.toString())

	// if the key type is known, set it explicitly
	// since older SSH server versions don't support key type negotiation,
	// and sometimes SSH clients hit the auth failure threshold before finding the right algorithm via brute force, see:
	// https://github.com/mwiede/jsch/issues/45#issuecomment-839727926
	if (keytype != null) {
		JSch.setConfig("PubkeyAcceptedAlgorithms", keytype)
	}

	// capture the JSch log
	val log = StringBuilder()
	JSch.setLogger(object : JschLogger {
		override fun isEnabled(level: Int) = true
		override fun log(level: Int, message: String?) {
			val levelstr = when (level) {
				JschLogger.DEBUG -> "DEBUG"
				JschLogger.INFO -> "INFO"
				JschLogger.WARN -> "WARN"
				JschLogger.ERROR -> "ERROR"
				JschLogger.FATAL -> "FATAL"
				else -> "???"
			}
			log.append("\t$levelstr: $message\n")
		}
	})

	// open the SSH connection
	val session = jsch.getSession(user, host, port)
	session.setConfig("PreferredAuthentications", "publickey")
	try {
		session.connect()
	} catch (t: Throwable) {
		System.err.println("""
			|Error connecting to SSH server. Troubleshooting tips:
			|   Make sure your username is correct: $user  (the username should not be quoted)
			|   Check that the private key is correct: $keypriv (exists? ${Files.exists(keypriv)})
			|      If your private key is not at the default location, set `dlab.key.private = /path` in gradle.properties.
			|   Check that the public key is correct: $keypub (exists? ${Files.exists(keypub)})
			|      If your public key is not at the default location, set `dlab.key.public = /path` in gradle.properties.
			|   Make sure the SSH client can find your key type before the auth failure limit (sometimes as low as 2).
			|      Check the SSH log (below) for details.
			|      If the correct key type isn't found before the auth failure limit, try setting `dlab.key.type = your-key=type` in gradle.properties.
			|      For older SSH keys, the correct type is often `ssh-rsa`.
			|   Check that the known hosts file is correct: $knownHosts (exists? ${Files.exists(knownHosts)})
			|      If your known_hosts file is not at the default location, set `dlab.knownHosts = /path` in gradle.properties.
			|   Make sure the SSH host is in your known_hosts file. Try connecting to $host via `ssh` in a terminal first.
			|   You'll need a Duke CS account to connect to the Duke CS SSH server.
			|SSH log:
			|$log
		""".trimMargin())
		throw t
	}
	try {
		val channel = session.openChannel("sftp") as ChannelSftp
		try {
			channel.connect()

			// at long last, we're connected. Do The Thing!
			return channel.block()

		} finally {
			channel.disconnect()
		}
	} finally {
		if (session.isConnected) {
			session.disconnect()
		}
	}
}

class SftpProgressLogger : SftpProgressMonitor {

	private var bytes: Long? = null
	private var progress: Long = 0
	private var startNs: Long? = null
	private var lastNs: Long? = null

	override fun init(op: Int, src: String?, dest: String?, max: Long) {
		println("Copying $max bytes from $src to $dest ...")
		bytes = max
		startNs = System.nanoTime()
		lastNs = startNs
	}

	override fun count(count: Long): Boolean {

		progress += count

		// don't log more than once per second
		val nowNs = System.nanoTime()
		val lastNs = lastNs
		val elapsedNs = if (lastNs != null) {
			nowNs - lastNs
		} else {
			0
		}
		if (elapsedNs >= 1_000_000_000) {

			this.lastNs = nowNs

			val bytes = bytes
			if (bytes != null) {
				val pct = 100f*progress.toFloat()/bytes.toFloat()
				println("   copied $progress bytes ${"%.1f".format(pct)} %")
			} else {
				println("   copied $progress bytes")
			}
		}

		return true // true to continue, false to cancel
	}

	override fun end() {
		val nowNs = System.nanoTime()
		val startNs = startNs
		if (startNs != null) {
			println("   done in ${"%.1f".format((nowNs - startNs).toFloat()/1e9f)} s")
		} else {
			println("   done")
		}
	}
}
