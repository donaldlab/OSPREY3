package osprey


object Jvm {

	const val javaLangVersion = 16

	const val packagePath = "edu/duke/cs/osprey"

	// add the module dependencies directly to the javac args
	// I don't think gradle has a good way to handle this yet?
	val moduleArgs = listOf(
		"--add-modules=jdk.incubator.foreign"
	)

	fun addModuleArgs(args: MutableList<String>?) {

		args ?: throw Error("args is null")

		// add the module args, if not already there
		args.addAll(moduleArgs.filter { it !in args })
	}
}
