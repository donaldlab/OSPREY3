package build

import edu.duke.cs.osprey.service.write
import org.jetbrains.dokka.CoreExtensions
import org.jetbrains.dokka.DokkaConfiguration.DokkaSourceSet
import org.jetbrains.dokka.base.DokkaBase
import org.jetbrains.dokka.base.resolvers.local.LocationProvider
import org.jetbrains.dokka.base.translators.documentables.DefaultDocumentableToPageTranslator
import org.jetbrains.dokka.links.DRI
import org.jetbrains.dokka.model.*
import org.jetbrains.dokka.model.doc.DocTag
import org.jetbrains.dokka.model.doc.DocumentationNode
import org.jetbrains.dokka.model.doc.Text
import org.jetbrains.dokka.pages.*
import org.jetbrains.dokka.plugability.*
import org.jetbrains.dokka.renderers.Renderer
import org.jetbrains.dokka.transformers.documentation.DocumentableToPageTranslator
import org.json.JSONArray
import org.json.JSONObject


// see: https://kotlin.github.io/dokka/1.6.0/developer_guide/introduction/

class OspreyPlugin : DokkaPlugin() {

	val dokkaBasePlugin by lazy { plugin<DokkaBase>() }

	val pager by extending {
		CoreExtensions.documentableToPageTranslator providing { ctx -> OspreyPager(ctx) } override dokkaBasePlugin.documentableToPageTranslator
	}

	val renderer by extending {
		CoreExtensions.renderer providing { ctx -> OspreyRenderer(ctx) } override dokkaBasePlugin.htmlRenderer
	}
}

private data class Config(
	val filename: String,
	val rootPackage: String
)

private fun DokkaContext.config(): Config {

	// read the plugin config
	val config = configuration.pluginsConfiguration
		.find { it.fqPluginName == OspreyPlugin::class.qualifiedName }
		?.values
		?: throw NoSuchElementException("need to send configuration JSON")
	val json = JSONObject(config)

	// read config values
	return Config(
		filename = json.getString("filename")!!,
		rootPackage = json.getString("package")!!
	)
}

private class OspreyPager(val ctx: DokkaContext): DocumentableToPageTranslator {

	class Page(
		override val name: String,
		val json: JSONObject
	) : RootPageNode() {

		override val children: List<PageNode>
			get() = emptyList()

		override fun modified(name: String, children: List<PageNode>): RootPageNode =
			Page(name, json)
	}

	private val config = ctx.config()
	private var locations: LocationProvider? = null

	override fun invoke(module: DModule): Page {

		// call the default pager to get URLs for all the elements
		val defaultPager = DefaultDocumentableToPageTranslator(ctx)
		val defaultPage = defaultPager.invoke(module)
		locations = ctx.plugin<DokkaBase>()
			.querySingle { locationProviderFactory }
			.getLocationProvider(defaultPage)

		// init the json with top-level elements
		val out = JSONObject()
		out.put("classlikes", JSONObject())
		out.put("funcs", JSONObject())
		out.put("props", JSONObject())

		// transform the documentation into JSON
		transform(module, out)

		// pass the json along the dokka pipeline as a page
		return Page("page for module ${module.dri}", out)
	}

	private fun urlKotlin(dri: DRI, sourceSets: Set<DokkaSourceSet>): String? {

		val displaySourceSets = sourceSets
			.map { it.toDisplaySourceSet() }
			.toSet()

		val locations = locations
			?: throw IllegalStateException("no locations available yet")

		return locations.resolve(dri, displaySourceSets)
	}

	private fun urlJava(dri: DRI): String {

		// assuming this DRI points to a java class or a class member, generate its javadoc URL
		// eg:      edu.duke.cs.osprey.structure/Molecule
		// becomes: edu/duke/cs/osprey/structure/Molecule.html

		val pname = dri.packageName
			?.replace('.', '/')
			?: throw Error("can't find URL to java element, no package: $dri")
		val cname = dri.classNames
			?: throw Error("can't find URL to java element, no class names: $dri")

		return "$pname.$cname.html"
	}

	private fun findClasslike(id: String, out: JSONObject): JSONObject? {
		val classlikes = out.getJSONObject("classlikes")
		return if (classlikes.has(id)) {
			classlikes.getJSONObject(id)
		} else {
			null
		}
	}

	private fun findClasslike(dri: DRI, out: JSONObject): JSONObject? =
		if (dri.classNames != null) {
			findClasslike(
				listOf(
					dri.packageName ?: "",
					dri.classNames
				).joinToString("/"),
				out
			)
		} else {
			null
		}

	private fun transform(doc: Documentable, out: JSONObject) =
		when (doc) {

			// ignore these
			is DTypeAlias,
			is DTypeParameter,
			is DParameter -> Unit

			// transform these
			is DModule -> transformModule(doc, out)
			is DPackage -> transformPackage(doc, out)
			is DEnumEntry -> transformEnumEntry(doc, out)
			is DFunction -> transformFunction(doc, out)
			is DClasslike -> transformClasslike(doc, out)
			is DProperty -> transformProperty(doc, out)

			// just in case
			else -> throw Error("don't know how to render Documentable: ${doc::class.qualifiedName}")
		}

	private fun transformModule(module: DModule, out: JSONObject) {

		// recurse
		for (doc in module.children) {
			transform(doc, out)
		}
	}

	private fun transformPackage(pack: DPackage, out: JSONObject) {

		// recurse
		for (doc in pack.children) {
			transform(doc, out)
		}
	}

	private fun transformClasslike(classlike: DClasslike, out: JSONObject) {

		// make a unique id for the class
		val id = listOf(
			classlike.dri.packageName
				?: "",
			classlike.dri.classNames
				?: throw NoSuchElementException("class ${classlike.name} has no class names... somehow")
		).joinToString("/")

		val outClass = JSONObject()
		out.getJSONObject("classlikes").put(id, outClass)

		// get the class name
		outClass.put("name", classlike.name)

		// get the KDoc, if any
		outClass.putIfSomething("kdoc", transformKdoc(classlike.documentation))

		// render the type info
		outClass.put("type", transformType(classlike.dri, classlike.sourceSets))

		// recurse
		for (doc in classlike.children) {
			transform(doc, out)
		}
	}

	private fun transformFunction(func: DFunction, out: JSONObject) {

		val outFunc = JSONObject()

		// get the function id and make a top-level object
		val baseId = listOf(
			func.dri.packageName ?: "",
			func.dri.classNames ?: "",
			func.name
		).joinToString("/")
		val signature =
			func.dri.callable?.signature() ?: ""
		val id = listOf(
			baseId,
			signature
		).joinToString("/")

		val outFuncs = out.getJSONObject("funcs")
		outFuncs.put(id, outFunc)

		// save the overload lookups
		outFuncs.getOrPutJSONArray(baseId)
			.put(signature)

		// add the function to the class, if any
		findClasslike(func.dri, out)
			?.getOrPutJSONArray("funcs")
			?.put(func.name)

		outFunc.put("name", func.name)

		// arguments
		val outArgs = JSONArray()
		for (param in func.parameters) {

			val outParam = JSONObject()

			outParam.put("name", param.name)
			outParam.put("type", transformTypeTree(param.type, param.sourceSets))

			outArgs.put(outParam)
		}
		outFunc.put("args", outArgs)

		// returns
		outFunc.put("returns", transformTypeTree(func.type, func.sourceSets))

		// receiver?
		val receiver = func.receiver
		if (receiver != null) {
			outFunc.put("receiver", transformTypeTree(receiver.type, receiver.sourceSets))
		}

		// link
		outFunc.put("url", "${config.rootPackage}/${urlKotlin(func.dri, func.sourceSets)}")

		// get the KDoc, if any
		outFunc.putIfSomething("kdoc", transformKdoc(func.documentation))
	}

	private fun transformProperty(prop: DProperty, out: JSONObject) {

		val outProp = JSONObject()

		// get the property id and make a top-level object
		val id = listOf(
			prop.dri.packageName ?: "",
			prop.dri.classNames ?: "",
			prop.name,
			prop.dri.callable?.signature() ?: ""
		).joinToString("/")
		out.getJSONObject("props").put(id, outProp)

		// add the property to the class, if any
		findClasslike(prop.dri, out)
			?.getOrPutJSONArray("props")
			?.put(prop.name)

		outProp.put("name", prop.name)
		outProp.put("type", transformTypeTree(prop.type, prop.sourceSets))

		// get the KDoc, if any
		outProp.putIfSomething("kdoc", transformKdoc(prop.documentation))
	}

	private fun transformEnumEntry(entry: DEnumEntry, out: JSONObject) {

		val outEntry = JSONObject()

		// look up the enclosing class
		val classId = listOf(
			entry.dri.packageName
				?: "",
			entry.dri.classNames
				?.let { names ->
					val suffix = ".${entry.name}"
					if (names.endsWith(suffix)) {
						names.substring(0, names.length - suffix.length)
					} else {
						throw IllegalStateException("malformed enum DRI: ${entry.name} not in ${entry.dri}")
					}
				}
				?: throw NoSuchElementException("enum ${entry.name} has no class names... somehow")
		).joinToString("/")
		val outClass = findClasslike(classId, out)
			?: throw NoSuchElementException("enum entry enclosing class not found: $classId")

		val outEntries = outClass.getOrPutJSONArray("entries")
		val outEntriesLut = outClass.getOrPutJSONObject("entriesLut")

		outEntry.put("name", entry.name)

		// get the KDoc, if any
		outEntry.putIfSomething("kdoc", transformKdoc(entry.documentation))

		outEntriesLut.put(entry.name, outEntries.length())
		outEntries.put(outEntry)
	}

	private fun transformKdoc(kdoc: SourceSetDependent<DocumentationNode>): String =
		kdoc.values.joinToString("\n") { node ->
			node.children.joinToString("\n") { child ->
				transformKdocTag(child.root)
			}
		}

	private fun transformKdocTag(tag: DocTag): String =
		when (tag) {

			// this appears to be the only tag with direct content
			is Text -> tag.body

			// everything else is a tree that may eventually have text children
			// NOTE: currently, no KDoc comments use tags, so I won't worry about trying
			// to render all the different tags into text correctly, for now
			else -> tag.children.joinToString("\n") {
				transformKdocTag(it)
			}
		}

	private fun transformTypeTree(type: Projection, sourceSets: Set<DokkaSourceSet>): JSONObject =
		when (type) {

			is GenericTypeConstructor -> transformType(type.dri, sourceSets).apply {

				// add generic parameters, if any
				if (type.projections.isNotEmpty()) {
					val outProjections = JSONArray()
					for (projection in type.projections) {
						outProjections.put(transformTypeTree(projection, sourceSets))
					}
					put("params", outProjections)
				}
			}

			is FunctionalTypeConstructor -> transformType(type.dri, sourceSets).apply {
				// TODO: do we need to bother with redenring functional types correctly?
				//  do they even appear in the Python API at all?
				put("functional", true)
			}

			is Nullable -> transformTypeTree(type.inner, sourceSets).apply {
				put("nullable", true)
			}

			is Variance<*> -> transformTypeTree(type.inner, sourceSets).apply {
				val variance = when (type) {
					is Covariance<*> -> "out"
					is Contravariance<*> -> "in"
					is Invariance<*> -> null
				}
				if (variance != null) {
					put("variance", variance)
				}
			}

			is TypeParameter -> JSONObject().apply {
				put("name", type.name)
			}

			is TypeAliased -> transformTypeTree(type.typeAlias, sourceSets)

			is Star -> JSONObject().apply {
				put("name", "*")
			}

			else -> throw Error("don't know how to render type: ${type::class.simpleName}")
		}

	private fun transformType(dri: DRI, sourceSets: Set<DokkaSourceSet>): JSONObject {

		val out = JSONObject()

		// get the dokka URL for this type, if any
		val url = if (dri.packageName?.startsWith(config.rootPackage) == true) {

			// one of ours, there definitely should be a url here
			val url = urlKotlin(dri, sourceSets)
				?: urlJava(dri)
			"${config.rootPackage}/$url"

		} else {

			// try to find the URL, dokka will still give us URLs for stdlib targets
			urlKotlin(dri, sourceSets)
		}

		if (url != null) {

			// have a url, render a short name
			out.put("url", url)
			out.put("name", dri.classNames?.split(".")?.last() ?: "(Unnamed)")

		} else {

			// no url, so use the full name
			out.put("name", "${dri.packageName}.${dri.classNames}")
		}

		return out
	}
}

private class OspreyRenderer(val ctx: DokkaContext) : Renderer {

	private val config = ctx.config()

	override fun render(root: RootPageNode) {

		// write out the file from the page JSON
		val page = root as OspreyPager.Page
		val file = ctx.configuration.outputDir.resolve(config.filename).toPath()
		page.json.toString(2).write(file)
	}
}


private fun JSONObject.getOrPutJSONArray(key: String): JSONArray =
	if (has(key)) {
		getJSONArray(key)
	} else {
		val new = JSONArray()
		put(key, new)
		new
	}

private fun JSONObject.getOrPutJSONObject(key: String): JSONObject =
	if (has(key)) {
		getJSONObject(key)
	} else {
		val new = JSONObject()
		put(key, new)
		new
	}

private fun JSONObject.putIfSomething(key: String, value: String?) {
	if (value != null && value.isNotBlank()) {
		put(key, value)
	}
}
