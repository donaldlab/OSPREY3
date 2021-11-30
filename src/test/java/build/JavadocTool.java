package build;


import com.sun.source.doctree.*;
import com.sun.source.tree.*;
import com.sun.source.util.*;
import edu.duke.cs.osprey.tools.FileTools;
import org.json.JSONArray;
import org.json.JSONObject;

import javax.lang.model.element.*;
import javax.lang.model.util.Elements;
import javax.tools.DiagnosticCollector;
import javax.tools.JavaFileObject;
import javax.tools.SimpleJavaFileObject;
import javax.tools.ToolProvider;
import java.io.IOException;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;


/**
 * Calls raw javac to get the full AST for source files,
 * because javadoc uses a truncated AST that doesn't have variable initializers anymore T_T
 */
public class JavadocTool {

	public static void main(String[] args) {

		// parse the args
		if (args.length != 2) {
			System.err.println("argument: <root package> <source folder or source file>");
			System.exit(1);
			return; // hahaha, stupid javac
		}

		var rootPackage = args[0];

		var path = Paths.get(args[1]);
		if (!Files.exists(path)) {
			System.err.println("source path does not exist: " + path);
			System.exit(1);
			return; // hahaha, stupid javac
		}

		var out = new JSONObject();

		int code;
		if (Files.isDirectory(path)) {
			code = runFolder(rootPackage, path, out);
		} else {
			code = runFile(rootPackage, path, out);
		}

		// write the JSON to stdout
		System.out.println(out.toString(2));
		System.exit(code);
	}

	public static int runFolder(String rootPackage, Path sourceFolder, JSONObject out) {

		int failed = 0;

		// walk the source folder
		List<Path> paths;
		try {
			paths = Files.walk(sourceFolder).toList();
		} catch (IOException e) {
			e.printStackTrace(System.err);
			return 2;
		}
		for (var path : paths) {

			// skip folders
			if (Files.isDirectory(path)) {
				continue;
			}

			try {
				handleFile(rootPackage, path.toString(), out);
			} catch (Exception ex) {
				System.err.println("can't process " + path);
				ex.printStackTrace(System.err);
				failed += 1;
			}
		}

		if (failed > 0) {
			System.err.println("failed to process " + failed + " files");
			return 3;
		}

		return 0;
	}

	public static int runFile(String rootPackage, Path sourceFile, JSONObject out) {

		try {
			handleFile(rootPackage, sourceFile.toString(), out);
		} catch (IOException ex) {
			System.err.println("can't process " + sourceFile);
			ex.printStackTrace(System.err);
			return 4;
		}

		return 0;
	}

	private static class Context {

		CompilationUnitTree tree;
		DocTrees trees;
		Elements elements;
		String rootPackage;

		TreePath path(Tree node) {
			return trees.getPath(tree, node);
		}

		Element element(Tree node) {
			return trees.getElement(path(node));
		}

		DocContext doc(Tree node) {

			var path = path(node);
			var docTree = trees.getDocCommentTree(path);
			if (docTree == null) {
				return null;
			}

			var ctx = new DocContext();
			ctx.tree = tree;
			ctx.trees = trees;
			ctx.elements = elements;
			ctx.rootPackage = rootPackage;
			ctx.docPath = path;
			ctx.docTree = docTree;
			return ctx;
		}

		String binaryName(ClassTree c) {
			return binaryName((TypeElement)element(c));
		}

		String binaryName(TypeElement elem) {
			return elements.getBinaryName(elem).toString();
		}
	}

	private static class DocContext extends Context {

		TreePath docPath;
		DocCommentTree docTree;

		DocTreePath path(DocTree node) {
			return DocTreePath.getPath(docPath, docTree, node);
		}

		Element element(DocTree node) {
			return trees.getElement(path(node));
		}
	}

	private static void handleFile(String rootPackage, String sourcePath, JSONObject out)
	throws IOException {

		// wrap the source in an object javac understands
		var sourceUrl = URI.create("string:///" + sourcePath);
		var sourceFile = new SimpleJavaFileObject(sourceUrl, JavaFileObject.Kind.SOURCE) {
			@Override
			public CharSequence getCharContent(boolean ignoreEncodingErrors) {
				return FileTools.readFile(sourcePath);
			}
		};

		// init the compiler
		var compiler = ToolProvider.getSystemJavaCompiler();
		var diagnostics = new DiagnosticCollector<JavaFileObject>();
		var task = (JavacTask)compiler.getTask(
			null,
			null,
			diagnostics,
			null,
			null,
			List.of(sourceFile)
		);

		// parse the source
		var tree = task.parse().iterator().next();

		// analyze the source to resolve the types
		task.analyze();

		// build the context
		var ctx = new Context();
		ctx.tree = tree;
		ctx.trees = DocTrees.instance(task);
		ctx.elements = task.getElements();
		ctx.rootPackage = rootPackage;

		// look for top-level classes
		for (var type : tree.getTypeDecls()) {
			switch (type.getKind()) {
				case CLASS -> handleClass(ctx, (ClassTree)type, out);
			}
		}
	}

	private static void handleClass(Context ctx, ClassTree c, JSONObject out) {

		var outClass = new JSONObject();
		out.put(ctx.binaryName(c), outClass);

		outClass.put("type", renderTypeTree(ctx, c));

		// get the class javadoc, if any
		var doc = ctx.doc(c);
		if (doc != null) {
			outClass.put("javadoc", renderJavadoc(doc));
		}

		// get the fields, methods, and inner classes
		var outFields = new JSONObject();
		outClass.put("fields", outFields);
		var outMethods = new JSONObject();
		outClass.put("methods", outMethods);
		for (var member : c.getMembers()) {
			switch (member.getKind()) {
				case VARIABLE -> handleField(ctx, (VariableTree)member, outFields);
				case METHOD -> handleMethod(ctx, (MethodTree)member, outMethods);
				case CLASS -> handleClass(ctx, (ClassTree)member, out);
			}
		}
	}

	private static void handleField(Context ctx, VariableTree field, JSONObject out) {

		var outField = new JSONObject();
		out.put(field.getName().toString(), outField);

		// add the type
		outField.put("type", renderTypeTree(ctx, field.getType()));

		// add the javadoc, if any
		var doc = ctx.doc(field);
		if (doc != null) {
			outField.put("javadoc", renderJavadoc(doc));
		}

		// add the intializer, if any
		var init = field.getInitializer();
		if (init != null) {
			outField.put("initializer", init.toString());
		}
	}

	private static void handleMethod(Context ctx, MethodTree method, JSONObject out) {

		var outMethod = new JSONObject();
		outMethod.put("name", method.getName().toString());

		// add the signature
		var buf = new StringBuilder();
		buf.append('(');
		for (var arg : method.getParameters()) {
			if (buf.length() > 1) {
				buf.append(',');
			}
			buf.append(arg.getType().toString());
		}
		buf.append(')');

		// add the return type, if any (constructors don't have return types)
		if (method.getReturnType() != null) {
			buf.append(method.getReturnType().toString());
			outMethod.put("returns", renderTypeTree(ctx, method.getReturnType()));
		}

		outMethod.put("signature", buf.toString());

		// add the javadoc, if any
		var doc = ctx.doc(method);
		if (doc != null) {
			outMethod.put("javadoc", renderJavadoc(doc));
		}

		// add method arguments
		var outArgs = new JSONArray();
		for (var arg : method.getParameters()) {
			var outArg = new JSONObject();

			outArg.put("name", arg.getName().toString());
			outArg.put("type", renderTypeTree(ctx, arg.getType()));

			outArgs.put(outArg);
		}
		outMethod.put("args", outArgs);

		// names aren't necessarily unique, so decorate with the signature
		out.put(method.getName().toString() + buf.toString(), outMethod);
	}

	private static Object renderTypeTree(Context ctx, Tree tree) {

		// get the element node for the tree node, if any
		var elem = ctx.element(tree);
		if (elem == null) {

			// no element, must be something simple like a primitive type
			// just render it directly to a string
			return tree.toString();
		}

		if (elem instanceof TypeElement typeElem) {

			var outType = renderType(ctx, typeElem);

			// look for type parameters
			if (tree instanceof ParameterizedTypeTree typeTree) {
				var params = typeTree.getTypeArguments();
				if (params != null && !params.isEmpty()) {
					outType.put("params", params.stream().map(p -> renderTypeTree(ctx, p)).collect(Collectors.toList()));
				}
			}

			return outType;

		} else if (elem instanceof TypeParameterElement) {
			// a type parameter, like Class<T>, just render to a string
			return elem.toString();
		}

		throw new Error("don't know how to render element for " + elem.getClass());
	}

	private static JSONObject renderType(Context ctx, TypeElement typeElem) {

		var outType = new JSONObject();
		outType.put("name", typeElem.getSimpleName().toString());

		// get the javadoc url for this type, if any
		if (typeElem.getQualifiedName().toString().startsWith(ctx.rootPackage + '.')) {

			outType.put("name", typeElem.getSimpleName().toString());
			outType.put("url", ctx.binaryName(typeElem)
				.replace('.', '/')
				.replace('$', '.')
				+ ".html");

		} else {

			// not one of our classes, no url so use the fully qualified name
			outType.put("name", typeElem.getQualifiedName().toString());
		}

		return outType;
	}

	// collect info to help consumers render javadoc tags
	private static Object renderJavadoc(DocContext doc) {

		var out = new JSONObject();

		BiConsumer<String,Object> append = (key, item) -> {
			if (!out.has(key)) {
				out.put(key, new JSONArray());
			}
			out.getJSONArray(key).put(item);
		};

		interface CustomTagRenderer {
			// no TriConsumer interface in the JVM runtime, so need a custom interface here to make the lambda
			void render(String name, DocTree node, List<? extends DocTree> contents);
		};
		CustomTagRenderer customTagRenderer = (name, node, contents) -> {
			switch (name) {

				case "cite" ->
					append.accept("citations", renderCitation(node, contents));

				case "warning", "warn" ->
					append.accept("warnings", renderWarning(node, contents));

				case "note" ->
					append.accept("notes", renderNote(node, contents));

				case "todo" ->
					append.accept("todos", renderTodo(node, contents));

				default ->
					throw new Error("unrecognized non-standard javadoc inline tag: " + name);
			}
		};

		doc.docTree.accept(new DocTreeScanner<Void,Void>() {

			@Override
			public Void visitParam(ParamTree node, Void ignored) {
				append.accept("params", renderParam(doc, node));
				return null;
			}

			@Override
			public Void visitLink(LinkTree node, Void ignored) {
				append.accept("links", renderLink(doc, node));
				return null;
			}

			// why, oh why wouldn't they let UnknownInlineTagTree and UnknownInlineBlockTree share a meaningful interface?!

			@Override
			public Void visitUnknownInlineTag(UnknownInlineTagTree node, Void ignored) {
				customTagRenderer.render(node.getTagName(), node, node.getContent());
				return null;
			}

			@Override
			public Void visitUnknownBlockTag(UnknownBlockTagTree node, Void ignored) {
				customTagRenderer.render(node.getTagName(), node, node.getContent());
				return null;
			}

		}, null);

		// get the full text, but just the body portion
		out.put("text", doc.docTree.getFullBody().stream()
			.map(n -> n.toString())
			.collect(Collectors.joining())
		);

		return out;
	}

	private static JSONObject renderParam(DocContext doc, ParamTree node) {

		var out = new JSONObject();

		// keep the raw string, for easy search-and-replace
		out.put("text", node.toString());

		// easy peasy
		out.put("name", node.getName().toString());
		out.put("description", node.getDescription().stream()
			.map(n -> n.toString())
			.collect(Collectors.joining())
		);

		return out;
	}

	private static JSONObject renderLink(DocContext doc, LinkTree node) {

		var out = new JSONObject();

		// keep the raw string, for easy search-and-replace
		out.put("text", node.toString());

		String label = null;

		// resolve the target type
		var elem = doc.element(node.getReference());
		if (elem == null) {

			throw new Error("unknown link target: " + node.getReference());

		} else if (elem instanceof TypeElement typeElem) {

			// class reference, just render the type
			out.put("type", renderType(doc, typeElem));

			// use the simple name, if we need a label
			if (label == null) {
				label = typeElem.getSimpleName().toString();
			}

		} else if (elem instanceof VariableElement varElem) {

			// variable reference, render the type of the containing element
			var container = varElem.getEnclosingElement();
			if (container instanceof TypeElement containerTypeElem) {
				out.put("type", renderType(doc, containerTypeElem));
			} else {
				throw new Error("unrecognized link target container: " + container.getClass().getSimpleName());
			}

			// and add the signature of the variable
			// but formatted for javadoc named anchors, ie:
			// fieldName
			out.put("signature", elem.getSimpleName().toString());

			// use the simple name, if we need a label
			if (label == null) {
				label = varElem.getSimpleName().toString();
			}

		} else if (elem instanceof ExecutableElement execElem) {

			// method/constructor reference, render the type of the containing element
			var container = execElem.getEnclosingElement();
			if (container instanceof TypeElement containerTypeElem) {
				out.put("type", renderType(doc, containerTypeElem));
			} else {
				throw new Error("unrecognized link target container: " + container.getClass().getSimpleName());
			}

			// and add the signature of the method
			// but formatted for javadoc named anchors, ie:
			// methodName(arg1,arg2)
			var buf = new StringBuilder();
			buf.append(elem.getSimpleName().toString());
			buf.append('(');
			int i = 0;
			for (var argElem : execElem.getParameters()) {
				if (i > 0) {
					buf.append(',');
				}
				i += 1;
				buf.append(argElem.asType().toString());
			}
			buf.append(')');
			out.put("signature", buf.toString());

			// use the simple name, if we need a label
			if (label == null) {
				label = execElem.getSimpleName().toString();
			}

		} else {
			throw new Error("unrecognized link target type: " + elem.getClass().getSimpleName());
		}

		// render the label
		if (!node.getLabel().isEmpty()) {
			label = node.getLabel().stream()
				.map(n -> n.toString())
				.collect(Collectors.joining());
		}
		out.put("label", label);

		return out;
	}

	private static JSONObject renderCitation(DocTree node, List<? extends DocTree> contents) {

		var out = new JSONObject();

		// keep the raw string, for easy search-and-replace
		out.put("text", node.toString());

		// preserve the line breaks in citations
		var lines = new JSONArray();
		for (var content : contents) {
			for (var line : content.toString().split("\n")) {
				lines.put(line.trim());
			}
		}
		out.put("lines", lines);

		return out;
	}

	private static JSONObject renderWarning(DocTree node, List<? extends DocTree> contents) {

		var out = new JSONObject();

		// keep the raw string, for easy search-and-replace
		out.put("text", node.toString());

		// just grab the first content node
		// TOOD: do we care about any other content nodes after that?
		var content = contents.stream().findFirst()
			.orElseThrow(() -> new Error("@warning has no content!"));
		out.put("content", content.toString());

		return out;
	}

	private static JSONObject renderNote(DocTree node, List<? extends DocTree> contents) {

		var out = new JSONObject();

		// keep the raw string, for easy search-and-replace
		out.put("text", node.toString());

		// just grab the first content node
		// TOOD: do we care about any other content nodes after that?
		var content = contents.stream().findFirst()
			.orElseThrow(() -> new Error("@note has no content!"));
		out.put("content", content.toString());

		return out;
	}

	private static JSONObject renderTodo(DocTree node, List<? extends DocTree> contents) {

		var out = new JSONObject();

		// keep the raw string, for easy search-and-replace
		out.put("text", node.toString());

		// just grab the first content node
		// TOOD: do we care about any other content nodes after that?
		var content = contents.stream().findFirst()
			.orElseThrow(() -> new Error("@todo has no content!"));
		out.put("content", content.toString());

		return out;
	}
}
