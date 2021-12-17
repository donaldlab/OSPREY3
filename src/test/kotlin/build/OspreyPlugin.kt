package build

import edu.duke.cs.osprey.service.write
import org.jetbrains.dokka.CoreExtensions
import org.jetbrains.dokka.DokkaConfiguration.DokkaSourceSet
import org.jetbrains.dokka.base.DokkaBase
import org.jetbrains.dokka.base.renderers.html.NavigationPage
import org.jetbrains.dokka.base.resolvers.local.LocationProvider
import org.jetbrains.dokka.links.DRI
import org.jetbrains.dokka.model.*
import org.jetbrains.dokka.model.doc.DocTag
import org.jetbrains.dokka.model.doc.DocumentationNode
import org.jetbrains.dokka.model.doc.Text
import org.jetbrains.dokka.pages.*
import org.jetbrains.dokka.plugability.*
import org.jetbrains.dokka.renderers.Renderer
import org.json.JSONArray
import org.json.JSONObject


// see: https://kotlin.github.io/dokka/1.6.0/developer_guide/introduction/

class OspreyPlugin : DokkaPlugin() {

	val dokkaBasePlugin by lazy { plugin<DokkaBase>() }

	val ospreyPlugin by extending {
		CoreExtensions.renderer providing { ctx -> OspreyRenderer(ctx) } override dokkaBasePlugin.htmlRenderer
	}
}

class OspreyRenderer(val dokkaCtx: DokkaContext) : Renderer {

	data class Context(
		val rootPackage: String,
		val locations: LocationProvider
	)

	override fun render(root: RootPageNode) {

		// read the plugin config
		val config = dokkaCtx.configuration.pluginsConfiguration
			.find { it.fqPluginName == OspreyPlugin::class.qualifiedName }
			?.values
			?: throw NoSuchElementException("need to send configuration JSON")
		val json = JSONObject(config)

		// read config values
		val filename = json.getString("filename")!!
		val rootPackage = json.getString("package")!!

		// build the output file path
		val file = dokkaCtx.configuration.outputDir.resolve(filename).toPath()

		// build the context
		val ctx = Context(
			rootPackage,
			dokkaCtx.plugin<DokkaBase>().querySingle { locationProviderFactory }.getLocationProvider(root)
		)

		// TEMP: apply the usual preprocessors
		//val preprocessors = dokkaCtx.plugin<DokkaBase>().query { htmlPreprocessors }
		//val newRoot = preprocessors.fold(root) { acc, t -> t(acc) }

		// render the documentation into JSON
		val out = JSONObject()
		render(root, ctx, out)

		// write out the JSON file
		out.toString(2).write(file)
	}

	private fun render(node: PageNode, ctx: Context, out: JSONObject) {

		// render this page, depending on the type
		when (node) {

			is ContentPage -> renderContent(node, ctx, out)

			is RendererSpecificRootPage,
			is RendererSpecificResourcePage,
			is NavigationPage -> Unit // ignore

			else -> throw Error("don't know how to render page ${node::class.qualifiedName}")
		}

		// recurse down the page tree
		for (child in node.children) {
			render(child, ctx, out)
		}
	}

	private fun renderContent(node: ContentPage, ctx: Context, out: JSONObject) =
		when (val doc = node.documentable) {

			// ignore these
			null,
			is DTypeAlias,
			is DPackage,
			is DModule,
			is DTypeParameter,
			is DParameter -> Unit

			// render these
			is DEnumEntry -> renderEnumEntry(doc, ctx, out)
			is DFunction -> renderFunction(doc, ctx, out)
			is DClasslike -> renderClasslike(doc, ctx, out)
			is DProperty -> renderProperty(doc, ctx, out)

			// just in case
			else -> throw Error("don't know how to render Documentable: ${doc::class.qualifiedName}")
		}

	private fun renderClasslike(classlike: DClasslike, ctx: Context, out: JSONObject) {

		// make a unique id for the class
		val id = listOf(
			classlike.dri.packageName
				?: "",
			classlike.dri.classNames
				?: throw NoSuchElementException("class ${classlike.name} has no class names... somehow")
		).joinToString("/")

		val outClass = JSONObject()
		out.put(id, outClass)

		// get the class name
		outClass.put("name", classlike.name)

		// get the KDoc, if any
		outClass.putIfSomething("kdoc", renderKdoc(classlike.documentation))

		// render the type info
		outClass.put("type", renderType(classlike, ctx))

		// TODO: properties, lut
		// TODO: functions, lut
	}

	private fun renderFunction(func: DFunction, ctx: Context, out: JSONObject) {

		// TODO
		//println("func ${func.dri}")
	}

	private fun renderProperty(prop: DProperty, ctx: Context, out: JSONObject) {

		val outProp = JSONObject()

		// look up the enclosing class, if any
		val classNames = prop.dri.classNames
		if (classNames != null) {

			val classId = listOf(
				prop.dri.packageName
					?: "",
				classNames
			).joinToString("/")
			val outClass = out.getJSONObject(classId)

			val outProps = outClass.getOrPutJSONArray("props")
			val outPropsLut = outClass.getOrPutJSONObject("propsLut")

			// add the property to the class
			outPropsLut.put(prop.name, outProps.length())
			outProps.put(outProp)

		} else {

			// otherwise, add it to the top-level object
			val id = listOf(
				prop.dri.packageName,
				"",
				prop.name
			).joinToString("/")
			out.put(id, outProp)
		}

		outProp.put("name", prop.name)
		outProp.put("type", renderTypeTree(prop.type, prop.sourceSets, ctx))
	}

	private fun renderEnumEntry(entry: DEnumEntry, ctx: Context, out: JSONObject) {

		val outEnum = JSONObject()

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
		val outClass = out.getJSONObject(classId)
			?: throw NoSuchElementException("enum entry enclosing class not found: $classId")

		val outEntries = outClass.getOrPutJSONArray("entries")
		val outEntriesLut = outClass.getOrPutJSONObject("entriesLut")

		val outEntry = JSONObject()

		outEntry.put("name", entry.name)

		// get the KDoc, if any
		outEntry.putIfSomething("kdoc", renderKdoc(entry.documentation))

		outEntriesLut.put(entry.name, outEntries.length())
		outEntries.put(outEntry)
	}

	private fun renderKdoc(kdoc: SourceSetDependent<DocumentationNode>): String =
		kdoc.values.joinToString("\n") { node ->
			node.children.joinToString("\n") { child ->
				renderKdocTag(child.root)
			}
		}

	private fun renderKdocTag(tag: DocTag): String =
		when (tag) {

			// this appears to be the only tag with direct content
			is Text -> tag.body

			// everything else is a tree that may eventually have text children
			// NOTE: currently, no KDoc comments use tags, so I won't worry about trying
			// to render all the different tags into text correctly, for now
			else -> tag.children.joinToString("\n") {
				renderKdocTag(it)
			}
		}

	private fun renderTypeTree(type: Projection, sourceSets: Set<DokkaSourceSet>, ctx: Context): JSONObject =
		when (type) {

			is GenericTypeConstructor -> renderType(type.dri, sourceSets, ctx).apply {

				// add generic parameters, if any
				if (type.projections.isNotEmpty()) {
					val outProjections = JSONArray()
					for (projection in type.projections) {
						outProjections.put(renderTypeTree(projection, sourceSets, ctx))
					}
					put("params", outProjections)
				}
			}

			is FunctionalTypeConstructor -> renderType(type.dri, sourceSets, ctx).apply {
				// TODO: do we need to bother with redenring functional types correctly?
				//  do they even appear in the Python API at all?
				put("functional", true)
			}

			is Nullable -> renderTypeTree(type.inner, sourceSets, ctx).apply {
				put("nullable", true)
			}

			is Variance<*> -> renderTypeTree(type.inner, sourceSets, ctx).apply {
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

			is TypeAliased -> renderTypeTree(type.typeAlias, sourceSets, ctx)

			else -> throw Error("don't know how to render type: ${type::class.simpleName}")
		}

	private fun renderType(classlike: DClasslike, ctx: Context): JSONObject =
		renderType(
			classlike.dri,
			classlike.sourceSets,
			ctx
		)

	private fun renderType(dri: DRI, sourceSets: Set<DokkaSourceSet>, ctx: Context): JSONObject {

		val out = JSONObject()

		// get the dokka URL for this type, if any
		if (dri.packageName?.startsWith(ctx.rootPackage) == true) {

			val displaySourceSets = sourceSets
				.map { it.toDisplaySourceSet() }
				.toSet()
			val url = ctx.locations.resolve(dri, displaySourceSets)
				?: throw NoSuchElementException("no URL found for $dri")
			out.put("url", "${ctx.rootPackage}/$url")

			// render a short name
			out.put("name", dri.classNames?.split(".")?.last() ?: "(Unnamed)")

		} else {

			// not one of our classes, no url so use the fully qualified name
			out.put("name", "${dri.packageName}.${dri.classNames}")
		}

		return out
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
