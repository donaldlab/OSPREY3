package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import cuchaz.kludge.window.FileDialog
import cuchaz.kludge.window.FilterList
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.gui.infoTip
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.tools.LZMA2
import edu.duke.cs.osprey.gui.compiler.*
import edu.duke.cs.osprey.gui.features.components.Icon
import edu.duke.cs.osprey.gui.features.components.icon
import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.forcefield.ForcefieldParams
import edu.duke.cs.osprey.gui.forcefield.amber.AmberForcefieldParams
import edu.duke.cs.osprey.gui.forcefield.eef1.EEF1ForcefieldParams
import edu.duke.cs.osprey.gui.io.UserSettings
import edu.duke.cs.osprey.gui.io.toBytes
import edu.duke.cs.osprey.gui.io.write
import edu.duke.cs.osprey.gui.prep.ConfSpace
import java.math.BigInteger
import java.nio.file.Path
import java.text.NumberFormat
import java.util.*


class CompileConfSpace(val confSpace: ConfSpace) : SlideFeature {

	override val id = FeatureId("compile.confspace")

	companion object {
		const val textExtension = "ccs"
		const val compressedExtension = "ccsx"
	}

	// make the compiler and configure the default settings
	private val compiler = ConfSpaceCompiler(confSpace).apply {
		forcefields.add(Forcefield.Amber96)
		forcefields.add(Forcefield.EEF1)
	}

	private val bigIntFormatter = NumberFormat.getIntegerInstance()
		.apply {
			isGroupingUsed = true
		}
	private fun BigInteger.format() =
		bigIntFormatter.format(this)

	private val winState = WindowState()
	private val pCompressed = Ref.of(true)

	private var progress: CompilerProgress? = null
	private var extraInfoBuf: Commands.TextBuffer = Commands.TextBuffer(1024)

	/** Cache the display info about the report, so we don't have to re-compute it every frame */
	private class ReportInfo(val report: ConfSpaceCompiler.Report) {

		val errorExtra: String? = report.error?.let { error ->

			// collect any extra information info a big buffer
			val msg = StringBuilder()

			report.error.extraInfo?.let { extra ->
				msg.append(extra)
				msg.append("\n")
			}

			error.causes().forEach { cause ->
				msg.append(cause.msgOrName())
				msg.append("\n")
			}

			msg
				.takeIf { it.isNotEmpty() }
				?.toString()
		}

		// TODO: save cause stack traces somewhere?
	}
	private var reportInfos = IdentityHashMap<ConfSpaceCompiler.Report,ReportInfo>()


	private fun resizeExtraInfoBufIfNeeded(msg: String): Commands.TextBuffer {

		if (!extraInfoBuf.hasRoomFor(msg)) {

			// make a bigger buffer
			extraInfoBuf = Commands.TextBuffer.of(msg)
		} else {

			// re-use the existing buffer
			extraInfoBuf.text = msg
		}

		return extraInfoBuf
	}

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Compile Conformation Space")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		winState.render(
			whenOpen = {

				// draw the window
				window("Compiler##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					tabBar("tabs") {
						tabItem("Conf Space") {
							guiConfSpace(imgui)
						}
						tabItem("Forcefields") {
							guiForcefields(imgui)
						}
						tabItem("Net Charges") {
							guiNetCharges(imgui)
						}
					}

					spacing()
					spacing()
					separator()
					spacing()
					spacing()

					// show either the compile button, or compilation progress
					val progress = progress
					if (progress == null) {

						// no compile running already,
						// show a button to start the compilation
						if (button("Compile")) {
							compile()
						}

					} else {

						val report = progress.report
						if (report == null) {

							// compilation running, show progress
							for ((taski, task) in progress.tasks.withIndex()) {

								if (taski > 0) {
									spacing()
								}

								progressBar(task.fraction, overlay = task.name)
							}

						} else {

							// compilation finished, show the report
							val reportInfo = reportInfos.getOrPut(report) { ReportInfo(report) }
							val compiled = report.compiled
							if (compiled != null) {
								guiSuccess(imgui, slidewin, reportInfo, compiled)
							} else {
								guiError(imgui, slidewin, reportInfo)
							}

							spacing()
							spacing()
							separator()
							spacing()
							spacing()

							// add a button to close the compilation report
							if (button("Finished")) {
								this@CompileConfSpace.progress = null
							}
						}
					}
				}
			}
		)
	}

	private fun guiConfSpace(imgui: Commands) = imgui.run {

		// show arbitrary stats about the conf space
		text("Conformation Space Info:")
		child("confSpaceInfo", 300f, 60f, true) {
			columns(2)

			text("Design Positions:")
			nextColumn()
			text("${confSpace.designPositionsByMol.values.sumBy { it.size }}")
			nextColumn()

			text("Sequences:")
			nextColumn()
			text(confSpace.positionConfSpaces.sequenceSpaceSize().format())
			nextColumn()

			text("Conformations:")
			nextColumn()
			text(confSpace.positionConfSpaces.confSpaceSize().format())
			nextColumn()

			columns(1)
		}
	}

	private var selectedForcefield: ForcefieldParams? = null

	private fun guiForcefields(imgui: Commands) = imgui.run {

		text("Forcefields:")
		child("forcefields", 200f, 100f, true) {
			for (ff in compiler.forcefields) {
				if (selectable(ff.forcefield.name, selectedForcefield === ff)) {
					selectedForcefield = ff
				}
			}
		}

		// if a forcefield is selected, show its settings
		selectedForcefield?.let { ff ->

			spacing()
			spacing()
			spacing()

			text("${ff.forcefield.name} settings:")
			indent(20f)

			when (ff) {

				is AmberForcefieldParams -> {

					inputDouble("Dielectric", Ref.of(ff::dielectric))
					sameLine()
					infoTip("""
						|The dielectric constant of the environment (aka its relative permittivity),
						|which influences electrostatic calculations
					""".trimMargin())

					checkbox("Distance-dependent dielectric", Ref.of(ff::distanceDependentDielectric))
					sameLine()
					infoTip("""
						|If checked, multiply the dielectric contstant by the atom pair distance (r)
						|for electrostatic interactions.
					""".trimMargin())

					inputDouble("van der Waals Scale", Ref.of(ff::vdwScale))
					sameLine()
					infoTip("""
						|Scaling to apply to the van der Waals calculations
					""".trimMargin())

					inputEnum("Charge Method", Ref.of(ff::chargeMethod))
					sameLine()
					infoTip("""
						|Method to generate partial charges for small molecules
					""".trimMargin())

					inputInt("Partial charge minimization steps", Ref.of(ff::sqmMinimizationSteps))
					sameLine()
					infoTip("""
						|Number of steps of minimization to perform during partial charge
						|calculation for small molecules. The default of 0 assumes the input
						|structures should not be perturbed while computing partial charges.
					""".trimMargin())
				}

				is EEF1ForcefieldParams -> {

					inputDouble("Scale", Ref.of(ff::scale))
					sameLine()
					infoTip("""
						|Scaling to apply to the solvent forcefield energy
					""".trimMargin())
				}

				else -> Unit
			}

			unindent(20f)
		}
	}

	private inner class NetChargeInfo(val mol: Molecule) {

		val label = mol.toString()
		// TODO: text storage

		// TODO: storage for each mutation

		val wildType: Ref<Int> = Ref.of(
			getter = { compiler.netCharges[mol]?.netCharge ?: 0 },
			setter = { compiler.netCharges[mol]?.netCharge = it }
		)

		init {
			// for the wild-type value, initialize from the molecule, or zero
			wildType.value = mol.netCharge ?: 0
		}
	}
	private val netChargeInfos =
		confSpace.mols
			.filter { (type, _) -> NetCharges.required(type) }
			.associate { (_, mol) -> mol to NetChargeInfo(mol) }

	private fun guiNetCharges(imgui: Commands) = imgui.run {

		for ((_, mol) in confSpace.mols) {
			val info = netChargeInfos[mol] ?: continue

			// breathe a little
			spacing()
			spacing()
			spacing()

			text(info.label)
			child(info.label, 300f, 100f, border = true) {

				inputInt("Wild Type", info.wildType)
				sameLine()
				infoTip("""
					|Input the formal net charge for this molecule before any mutations.
				""".trimMargin())
			}
		}
	}

	private fun compile() {
		progress = compiler.compile()
	}

	private fun guiSuccess(imgui: Commands, slidewin: SlideCommands, reportInfo: ReportInfo, compiled: CompiledConfSpace) = imgui.run {

		// TODO: formatting? make it pretty?

		text("Conformation space compiled successfully!")

		// TODO: render warning texts into the reportInfo
		// show compiler warnings
		val warnings = reportInfo.report.warnings
		if (warnings.isNotEmpty()) {
			text("Compiler warnings: ${warnings.size}")
			indent(20f)
			for (warning in warnings) {

				// show the warning
				icon(slidewin, Icon.SignWarning)
				sameLine()
				text(warning.msg)

				// show a button to show more info, if needed
				warning.extraInfo?.let { extra ->
					// TODO: show this somehow
				}
			}
			unindent(20f)
		}

		// show the save button
		if (button("Save")) {

			val currentExtension =
				if (pCompressed.value) {
					compressedExtension
				} else {
					textExtension
				}

			val filterList = FilterList(listOf(currentExtension))
			FileDialog.saveFile(filterList, UserSettings.openSaveDir)?.let { chosenPath ->
				UserSettings.openSaveDir = chosenPath.parent

				// decide if we should add an extension or not
				val filename = chosenPath.fileName.toString()
				val savePath = if (!filename.contains(".")) {

					// the user didn't choose an extension, so add one automatically
					chosenPath.parent.resolve("$filename.$currentExtension")

				} else {

					// otherwise, use exactly the filename the user gave us
					chosenPath
				}

				compiled.save(pCompressed.value, savePath)
			}
		}

		// show compression options
		sameLine()
		checkbox("Compress", pCompressed)
		sameLine()
		infoTip("""
			|Compiled conformation spaces can take up quite a bit of space.
			|Using compression will make the files significantly smaller,
			|and speed up transfers if you want to copy the compiled
			|conformation space files anywhere or send them over a network.
		""".trimMargin())
	}

	private fun CompiledConfSpace.save(compress: Boolean, path: Path) {

		// actually write out the compiled conf space to a file
		var bytes = toBytes()
		if (compress) {
			bytes = LZMA2.compress(bytes)
		}
		bytes.write(path)
	}

	private fun guiError(imgui: Commands, slidewin: SlideCommands, reportInfo: ReportInfo) = imgui.run gui@{

		// TODO: formatting? make it pretty?

		text("Conformation space failed to compile.")

		val error = reportInfo.report.error ?: run {
			// this shouldn't happen, but just in case ...
			text("but no error information was available")
			return@gui
		}

		// show the simple error message
		icon(slidewin, Icon.SignError)
		sameLine()
		text(error.msgOrName())

		// show the extra info, if needed
		reportInfo.errorExtra?.let { extra ->
			inputTextMultiline(
				"",
				resizeExtraInfoBufIfNeeded(extra),
				600f, 400f,
				IntFlags.of(Commands.InputTextFlags.ReadOnly)
			)
		}
	}
}
