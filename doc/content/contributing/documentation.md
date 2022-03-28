+++
menuTitle = "Documentation"
title = "OSPREY's Documentation System"
weight = 2
+++

## Overview

Osprey is mostly implemented in Java and Kotlin, but exposes a Python API as the main user interface for the
Server component. For the most part, this Python API is a very thin wrapper around the Java and Kotlin code
and much of it is dynamically generated using [JPype][jpype]. The goal of this arrangement is use the JVM and JVM
languages to enable rapid development of high-performance code, but provide end users with the ease and flexibility
of Python scripting. As much as possible, OSPREY's Server component tries to look and feel like a real Python module,
even though it's really not.

OSPREY's documentation exists in two main places:
 * Documents written in [Markdown][markdown]
 * Comments in code source files

[jpype]: https://jpype.readthedocs.io/en/latest/
[markdown]: https://en.wikipedia.org/wiki/Markdown

All of these files are processed using automated tools to produce
the website you're reading right now.
In the case of the Markdown documents, OSPREY uses [Hugo][hugo], a static website generator,
to convert the Markdown into HTML/CSS/Javascript for the website.

[hugo]: https://gohugo.io

For the comments in code source files, language-specific tools are used to convert the code comments
into either HTML/CSS/Javascript files directly, or into Markdown files that get sent to Hugo for conversion into
HTML/CSS/Javascript in subsequent processing steps.


### Java comments

Java comments are processed with the standard [JavaDoc][javadoc] tool directly into HTML/CSS/Javascript.

[javadoc]: https://docs.oracle.com/en/java/javase/17/javadoc/javadoc.html


### Kotlin comments

Kotlin comments are processed with the standard [Dokka][dokka] tool directly into HTML/CSS/Javascript.

[dokka]: https://kotlinlang.org/docs/kotlin-doc.html


### Python comments

Python, being an ecosystem with few real standards, has many options to choose from when it comes to
documentation. OSPREY uses a tool called [PyDoc Markdown][pydoc-markdown] to convert Python docstrings
into Markdown documents. These Markdown documents are then passed on to Hugo for conversion to HTML/CSS/Javascript.

[pydoc-markdown]: https://niklasrosenstein.github.io/pydoc-markdown/

In addition, OSPREY has customized the docstring processing in PyDoc Markdown to allow referencing documentation
and types from the Java and Kotlin code. Since OSPREY's Python API is often just a thin veneer over the Java
and Kotlin code, referencing other code types and documentation directly in the Python docstrings makes it easier to keep
the Python documentation consistent and synchronized with the underlying code in Java and Kotlin.

On the Python side, the customizations to PyDoc Markdown live in the `buildSrc/src/main/python` folder in the
`osprey_pydoc` module. The Python code there reads [JSON][json] metadata about the Java and Kotlin code to build
the documentation files from the docstrings in the Python source code.

[json]: https://www.json.org/json-en.html

On the Java side, OSPREY calls the Java compiler (in `buildSrc/src/main/java/osprey/JavadocTool.java`) to digest the Java
code into the JSON metadata.

On the Kotlin side, OSPREY calls the standard Dokka tool with a custom plugin
(in `buildSrc/src/main/kotlin/osprey/OspreyDokkaPlugin.kt`) to digest the Kotlin code into JSON metadata.


## Writing documentation

### in Markdown documents

OSPREY uses Hugo, with the [Learn theme][hugo-theme-learn], to render the Markdown documents into HTML/CSS/Javascript
for the website.

[hugo-theme-learn]: https://learn.netlify.app/en/

All the Markdown documents are stored in `doc/content` and organized into folders by category.
Each folder has a `_index.md` file in it which serves as the default (or "index") file for that folder.

Each Markdown document has a bit of non-Markdown code at the top, called [Front Matter][hugo-front-matter].
OSPREY uses the [Toml][toml] convention for Hugo front matter.
The front matter configures website options like:
 * which Entries appear in the Table of Contents menu on the left
 * what text shows on the link in the menu
 * what text shows as the title of the document
 * what order the documents should appear in the menu

[hugo-front-matter]: https://gohugo.io/content-management/front-matter/
[toml]: https://toml.io/en/

OSPREY's Hugo installation is configured (in `doc/config.toml`) to use the [Goldmark][goldmark] Markdown processor.
Goldmark implements the standardized [CommonMark][commonmark] dialect of Markdown, and also allows extensions
to the language. We've confiured Goldmark to also allow raw HTML to be rendered from Markdown documents.

[goldmark]: https://github.com/yuin/goldmark
[commonmark]: https://commonmark.org/

For more information about Hugo/Learn features useful for writing Markdown documents, see
[Hugo Learn Theme Documentation, Content Chapter](https://learn.netlify.app/en/cont/).


#### Hugo server mode

To help shorten the write-read-revise loop, Huge offers a server mode that watches for changes in the `doc` folder,
re-renders documents, and refreshes the browser page as needed. You can turn on server mode with the Gradle task
`hugoServer`. After starting, Hugo will print the URL of your temporary website. It might be something like
[http://localhost:1313/donaldlab/software/osprey/docs/](http://localhost:1313/donaldlab/software/osprey/docs/).
Point your browser to that URL to see quick feedback on your changes to Markdown files while you're editing them.
The Hugo server will keep running until you tell it to stop. You can stop the Hugo server by [stopping the running
Gradle task][gradle-stop-task].

[gradle-stop-task]: {{< ref "/contributing/gradle#stop-task" >}}


#### Some useful Markdown extensions available in OSPREY

Hugo provides [Shortcodes][hugo-shortcodes] to extend the Markdown syntax to implement useful features for documentation.

[hugo-shortcodes]: https://gohugo.io/content-management/shortcodes/

Here are some useful shortcodes for OSPREY:

#### Links to other Markdown documents

To link to another Markdown document, use the `ref` shortcode, with the path to the document relative to `/doc/content`
in the source tree.
```markdown
{{</* ref "/path/to/doc" */>}}
```
You can even link directly to sections in a document with a "hash":
```markdown
{{</* ref "/path/to/doc#hash" */>}}
```


#### Notice Boxes

You can draw attention to certain sections of the documentation using the `notice` shortcode.
For example:
```markdown
{{%/* notice tip */%}}
Try out this one weird tip that doctors hate!
{{%/* /notice */%}}
```
will render as:
{{% notice tip %}}
Try out this one weird tip that doctors hate!
{{% /notice %}}

Different kinds of notices are available too, like Note, Info, and Warning.
See the [Hugo Theme Learn Docuemntation][hugo-theme-learn-notice] for more information.

[hugo-theme-learn-notice]: https://learn.netlify.app/en/shortcodes/notice/


### in Java code

Documentation in Java code is written in the [JavaDoc][javadoc-lang] language.
And to be super confusing, the JavaDoc language shares the same name as the JavaDoc tool,
but the two are separate concepts.

[javadoc-lang]: https://www.oracle.com/technical-resources/articles/java/javadoc-tool.html

OSPREY has made a few small extensions to the JavaDoc language to allow more control over the appearance
of the eventual output in the Python code documentation on the website.

We've added a few new directives to the JavaDoc language:

 * `@cite` helps reference literature citations in the code documentation. For example:
    ```java
    /**
     * Here are some non-standard javadoc tags
     * {@cite This is a citation.
     * It can have multiple lines.
     * Sometimes even three.}
     */
   ```
   will render in the Python documentation as:
   {{% notice info %}}
   This is a citation.\
   It can have multiple lines.\
   Sometimes even three.
   {{% /notice %}}
 * `@warning` shows a warning to developers. For example:
   ```java
   /**
    * {@warning Be careful not to do the bad thing!}
    */
   ```
   will render in the Python documentation as:
   {{% notice warning %}}
   Be careful not to do the bad thing!
   {{% /notice %}}
 * `@note` shows a notice to developers. For example:
   ```java
   /**
    * {@note Notice this thing.}
    */
   ```
   will render in the Python documentation as:
   {{% notice note %}}
   Notice this thing.
   {{% /notice %}}
 * `@todo` tells developers that some feature has yet to be finished. For example:
   ```java
   /**
    * {@todo We haven't finished this thing yet.}
    */
   ```
   will render in the Python documentation as:
   {{% notice tip %}}
   **TODO:** We haven't finished this thing yet.
   {{% /notice %}}


### in Kotlin code

Documentation comments in Kotlin are written in the [KDoc][kdoc] language.

[kdoc]: https://kotlinlang.org/docs/kotlin-doc.html

KDoc is very similar to JavaDoc. However, the OSPREY extensions to JavaDoc
are currently not implemented for KDoc.


### in Python code

Code documentation in Python files takes the form of [docstrings][docstrings].
The docstrings are processed by PyDoc Markdown using the [PydocmdProcessor][pydoc-processor], which supports a
simple format. For example, for a function docstring:
```python
def my_func(arg1):
    """
    This is the description of my_func.
    
    # Arguments
    arg1 (int): The first argument.
    kwargs (dict): Keyword arguments.
    
    # Raises
    RuntimeError: If something bad happens.
    ValueError: If an invalid argument is specified.
    
    # Returns
    A value.
    """
```

[docstrings]: https://peps.python.org/pep-0257/
[pydoc-processor]: https://niklasrosenstein.github.io/pydoc-markdown/api/pydoc_markdown/processors/#class-pydocmdprocessor

In OSPREY, we have extensively extended PyDoc Markdown to reference information from the
Java and Kotlin documentation systems. The extensions take the form of a preprocessor macro system that expands
occurrences of `${cmd}` into PyDoc Markdown-compatible syntax, for various values of `cmd`.

At the current time of writing, over 20 different macro commands are supported.
The authoritative source for which commands are implemented is the Python source code itself,
at `buildSrc/src/main/python/osprey_pydoc/__init__.py`, in the `OspreyProcessor.__init__()` function.
The following is a brief description of each command, its arguments, and its effects.

#### `class_javadoc(path)`
Reproduces the JavaDoc content for the given [Java class identified by `path`](#classes-java).

**Examples:**
```
${class_javadoc(.parallelism.Parallelism)}
```
   
#### `class_kdoc(id)`
Reproduces the KDoc content for the given [Kotlin class identified by `id`](#classes-kotlin).

**Examples:**
```
${class_kdoc(.gui.forcefield/Forcefield.Amber96)}
```

#### `type_java(path)`
Renders the name of the [Java class identified by `path`](#classes-java) and links to the relevant JavaDoc page.

**Examples:**
```
${type_java(.structure.Molecule)}
```

#### `type_jre(path)`
Renders the name of the Java class identified by `path` in the Java runtime environment (JRE),
and links to the relevant JavaDoc page on the Oracle documentation website.
Unlike class paths in other macros, `path` here must be the fully-qualified name of the JRE class.

**Examples:**
```
${type_jre(java.util.ArrayList)}
```

#### `type_kotlin(id)`
Renders the name of the [Kotlin class identified by `id`](#classes-kotlin) and links to the relevant KDoc page.

**Examples:**
```
${type_kotlin(.gui.motions/DihedralAngle)}
```

#### `field_javadoc(path)`
Renders the JavaDoc content for the [Java field identified by `path`](#fields-java).

**Examples:**
```
${field_javadoc(.gmec.Comets$Builder#logFile)}
```

#### `arg_field_javadoc(arg_name, field_path, type=type)`
Renders an argument entry for a Python function that draws information from a Java class field.
Useful for Python API functions that wrap code written in the [Java "builder" class pattern][java-builder].
`arg_name` refers to the name of the argument in the Python function.
`field_path` is the [Java field](#fields-java).
`type` is an optional named argument to override the Java type with something more Pythonic, if needed.

[java-builder]: https://www.geeksforgeeks.org/builder-pattern-in-java/

**Examples:**
```
# Arguments
${arg_field_javadoc(shellDist, .confspace.SimpleConfSpace$Builder#shellDist)}
${arg_field_javadoc(maxNumNodes, .astar.conf.ConfAStarTree$Builder#maxNumNodes, type=int)}
```

#### `arg_fields_javadoc(path, [arg_name, field_name, type=type])`
Renders multiple arguments entries for a Python function where the information comes from a list of Java class fields.
`path` is the [Java class](#classes-java).
The second argument is an array of argument and field values:
`arg_name` names the argument in the Python function.
`field_name` names the field from the Java class identified by `path`. Note this is not a fully-qualified name for the field.
It is merely the name of the field in the Java class.
`type` is an optional named argument to override the Java type of the field with something more Pythonic, if needed.

**Examples:**
```
# Arguments
${args_fields_javadoc(.gmec.SimpleGMECFinder$Builder,
   [printIntermediateConfs, printIntermediateConfsToConsole],
   [useExternalMemory],
   [confDBFile, confDB, type=str]
)}
```

#### `returns_method_java(path)`
Renders the return type of the [Java method identified by `path`](#methods-java).

**Examples:**
```
# Returns
${returns_method_java(.pruning.SimpleDEE#read)}

# Returns
${returns_method_java(.structure.PDBIO#readFile(String)Molecule)}
```

#### `returns_func_kotlin(path)`
Renders the return type of the [Kotlin function identified by `id`](#functions-kotlin).

**Examples:**
```
# Returns
${returns_func_kotlin(.gui.forcefield.amber//findTypes)}

# Returns
${returns_func_kotlin(.gui.motions/DihedralAngle.ConfDescription.Companion/makeFromLibrary)}
```

#### `method_javadoc(path)`
Renders the JavaDoc of the [Java method identified by `path`](#methods-java).

**Examples:**
```
${method_javadoc(.astar.conf.ConfAStarTree$Builder#setTraditional)}
```

#### `func_kdoc()`
Renders the KDoc of the [Kotlin function identified by `id`](#functions-kotlin).

**Examples:**
```
${func_kdoc(.gui.forcefield.amber//inferProtonation)}
${func_kdoc(.gui.prep/Proteins/makeDesignPosition)}
```

#### `arg_java(python_name, method_path, java_name, type=type)`
Renders an argument entry for a Python function that draws information from a Java method argument.
`python_name` refers to the name of the argument in the Python function.
`method_path` is the [Java method](#methods-java).
`java_name` is an optional positional argument to give the name of the argument in the Java method,
if it differs from the name of the Python argument.
`type` is an optional named argument to override the Java type with something more Pythonic, if needed.


**Examples:**
```
# Arguments
${arg_java(kstar, .kstar.SequenceAnalyzer#<init>(KStar))} a configured instance of KStar
${arg_java(rcs, .astar.conf.RCs#<init>(RCs,PruningMatrix), other)}
```

#### `args_java(path, [python_name, java_name, type=type])`
Renders multiple arguments entries for a Python function where the information comes from a list of Java method arguments.
`path` is the [Java method](#methods-java).
The second argument is an array of function/method argument values:
`python_name` names the argument in the Python function.
`java_name` names the argument in the Java method identified by `path`.
`type` is an optional named argument to override the Java type with something more Pythonic, if needed.

**Examples:**
```
# Arguments
${args_java(.lute.LUTEPfunc#<init>,
    [luteEcalc, ecalc],
    [astar]
)}
${args_java(.sofea.Sofea$StateConfig#<init>,
    [emat],
    [confEcalc],
    [confdbPath, confDBFile, type=str]
)}
```

#### `arg_kotlin(python_name, function_id, kotlin_name, type=type)`
Renders an argument entry for a Python function that draws information from a Java method argument.
`python_name` refers to the name of the argument in the Python function.
`function_id` is the [Kotlin function](#functions-java).
`kotlin_name` is an optional positional argument to give the name of the argument in the Kotlin function,
if it differs from the name of the Python argument.
`type` is an optional named argument to override the Kotlin type with something more Pythonic, if needed.

**Examples:**
```
# Arguments
${arg_kotlin(mol, .gui.prep/DuplicateAtoms/DuplicateAtoms)}:
The molecule for which to search for duplicate atoms
${arg_kotlin(confSpace, .gui.compiler/ConfSpaceCompiler/ConfSpaceCompiler)}: the conformation space to compile
```

#### `args_kotlin(id, [python_name, kotlin_name, type=type])`
Renders multiple arguments entries for a Python function where the information comes from a list of Kotlin function arguments.
`id` is the [Kotlin function](#functions-kotlin).
The second argument is an array of function/method argument values:
`python_name` names the argument in the Python function.
`kotlin_name` names the argument in the Kotlin function identified by `id`.
`type` is an optional named argument to override the Java type with something more Pythonic, if needed.


**Examples:**
```
${args_kotlin(.molscope.molecule/Polymer/findResidueOrThrow/#kotlin.String#kotlin.String,
    [chainId],
    [resId]
)}
```

#### `receiver_kotlin(python_name, id)`
Renders the receiver of a [Kotlin function idenfitied by `id`](#functions-kotlin) as an argument to the Python function.
`python_name` names the argument in the Python function.

**Examples:**
```
# Arguments
${receiver_kotlin(mol, .gui.forcefield.amber//findTypes)}:
The molecule (or multiple molecules in a single Molecule instance) for which to determine molecule types

# Arguments
${receiver_kotlin(mol, .gui.forcefield.amber//deprotonate/.molscope.molecule.Molecule#)}:
The molecule from which to remove atoms
```

#### `arg_javadoc(path, arg_name)`
Renders the JavaDoc of the argument named `arg_name` in a [Java method identified by `path`](#methods-java).

**Examples:**
```
# Arguments
tempDir `str`: ${arg_javadoc(.externalMemory.ExternalMemory#setTempDir(String)void, dir)}
```

#### `enum_java(path)`
Renders the JavaDoc and all of the enum values for the [Java enum class identified by `path`](#classes-java).

**Examples:**
```
Precision = osprey.c.gpu.Structs.Precision
'''
${enum_java(.gpu.Structs$Precision)}
'''
```

#### `enum_kotlin(id)`
Renders the KDoc and all of the enum values for the [Kotlin enum class identified by `id`](#classes-kotlin).

**Examples:**
```
Hybridization = osprey.c.gui.forcefield.amber.Hybridization
'''
${enum_kotlin(.gui.forcefield.amber/Hybridization)}
'''
```

#### `default(arg_nane, value)`
Shows the default value of a Python function argument `arg_name` as the literal `value`.

**Examples:**
```
# Arguments
tempDir `str`: The temporary directory
${default(tempSubdir, <automatically generated>)}
```

#### `prop_kotlin(id)`
Renders the name, type, and KDoc of the [Kotlin property identified by `id`](#properties-kotlin).

**Examples:**
```
confLibs = osprey.c.gui.features.components.ConfLibs.INSTANCE.getInfos()
'''
${prop_kotlin(.gui.features.components/ConfLibs/infos)}
'''
```

#### `link_func_kotlin(id)`
Renders a link to the KDoc page for the` [Kotlin function identified by `id`](#functions-kotlin).

**Examples:**
```
${link_func_kotlin(.gui.prep/ConfSpace/getConformations/#edu.duke.cs.osprey.gui.prep.DesignPosition#kotlin.String)}
```
 

### Referring to Java entities with paths

#### Classes {#classes-java}

When a Java class is identified by `path`, the `path` should be the "binary name" of the class, 
which is similar to the [fully-qualified name of the class][java-fqn], except inner classes are delimited by
`$` instead of `.`. For example, the "binary name" of the JRE class `ArrayList` is `java.util.ArrayList`.
The "binary name" of the JRE class `Map.Entry` is `java.util.Map$Entry`.
As a shortcut, all OSPREY classes can be referenced with a relatively-qualified class name,
where the common prefix `edu.duke.cs.osprey` has been omitted and the remaining path begins with `.`.

**Examples:**
```
edu.duke.cs.osprey.package.Classname
.package.Classname
.package.Classname$Inner1
.package.Classname$Inner1$Inner2
```

[java-fqn]: https://docs.oracle.com/javase/specs/jls/se11/html/jls-6.html#jls-6.7


#### Fields {#fields-java}

When a Java field is identified by `path`, `path` follows the same rules as for paths that identify Java classes,
except the suffix `#fieldName` is appended to identify the field within the referenced class, where `fieldName`
names the field in the class.

**Examples:**
```
edu.duke.cs.osprey.package.Classname#field
.package.Classname#field
.package.Classname$Inner1#field
```


#### Methods {#methods-java}

When a Java method is identified by `path`, `path` follows the same rules as for paths that identify Java classes,
except the suffix `#methodName` is appended to identify the method within the referenced class, where `methodName`
names the method in the class.

Additionally, to disambiguate overloaded methods, you may need to append the method signature to the path,
e.g. `#methodName(arg1Type,arg2Type)returnType`.

It can be tricky to format the method signatures exactly right to resolve any overloads. OSPREY's extensions to
PyDoc Markdown attempt to help with this by offering suggestions when multiple overloads are detected.
If you start by specifying the method name without a method signature, OSPREY will throw a Python error with a list
of all of the possible method overloads and their signatures. Simply choose a method signature from one of the
overloads and try again.

For Java class constructors, the name `<init>` should be used.

**Examples:**
```
.pruning.SimpleDEE#read
.parallelism.Parallelism$Builder#build
.structure.PDBIO#readFile(String)Molecule
.paste.Paste#<init>
```


### Referring to Kotlin entities with ids

#### Classes {#classes-kotlin}

When a Kotlin class is identified by `id`, the `id` is the unique identifier used by Dokka for the class.
For Kotlin classes, `id` takes the form `$package/$classname` where:

 * `$package` is the package of the class.
 * `$classname` is the name of the class within the package.

As a shortcut, OSPREY classes can be referenced with a relatively-qualified package where the common
prefix `edu.duke.cs.osprey` has been omitted, and the rest of the package name begins with `.`.

**Examples:**
```
edu.duke.cs.osprey.gui.motions/DihedralAngle
.gui.motions/DihedralAngle
.gui.motions/DihedralAngle.LibrarySettings
```


#### Properties {#properties-kotlin}

When a Kotlin property is identified by `id`, the `id` is the unique identifier used by Dokka for the property.
For Kotlin properties, `id` takes the form `$package/$classname/$propertyName` where:

 * `$package` is the package of the class.
 * `$classname` is the name of the class within the package.
 * `$propertyName` is the name of the property within the class. 

As a shortcut, OSPREY classes can be referenced with a relatively-qualified package where the common
prefix `edu.duke.cs.osprey` has been omitted, and the rest of the package name begins with `.`.

**Examples:**
```
.gui.features.components/ConfLibs/infos
```


#### Functions {#functions-kotlin}

When a Kotlin function is identified by `id`, the `id` is the unique identifier used by Dokka for the function.
For Kotlin functions, `id` takes the form `$package/$className/$functionName/$signature` where:

 * `$package` is the package of the function, which can be relative or absolute (see [Kotlin classes](#classes-kotlin)).
 * `$className` is the name of the class within the package, if any. If the function is a free function and not a
    member of any class, the `$classname` is empty.
 * `$functionName` is the name of the function itself.
    For Kotlin class constructors, the name `<init>` should be used.
 * `$signature` is optional, but is necessary to disambiguate overloaded functions. It takes the form of:
   `#arg1Type#arg2Type`. If `$signature` is omitted, the preceding slash should also be omitted.

**Examples:**
```
.gui.prep/ConfSpace/getConformations/#edu.duke.cs.osprey.gui.prep.DesignPosition#kotlin.String
```


## Building the documentation

### Prerequisite tools

Before you can build the documentation, you'll need to install some tools:

[Hugo](https://gohugo.io):
 * [Download releases](https://github.com/gohugoio/hugo/releases)\
   Choose the "extended" version for your operating system. The extended features are needed by PyDoc Markdown.

[PyDoc Markdown](https://github.com/NiklasRosenstein/pydoc-markdown):
 * [Installation Instructions](https://niklasrosenstein.github.io/pydoc-markdown/#installation)\
   Make sure to install the "old-style" version with the YAML configuration file (before v4.6).
   Apparently PyDoc Markdown has been completely rewritten recently to use a new backend called Novella,
   but OSPREY hasn't yet been updated to use it. And why huge updates like this only seem to merit a minor
   version bump is beyond me... Open Source is fun like that I guess.


### Gradle tasks

To build the documentation, there are a few different tasks you'll use:

 * **`buildDocsRelease`**:\
   This task generates all the documentation that's specific to the current version of OSPREY,
   so it can be archived and shown on the website next to all the other versions of OSPREY.
   The task generates the documentation from the code comments, runs Hugo on all the Markdown documents
   in `/doc/content/documentation/main`, and packages the results into a release archive at
   `/build/releases/osprey-docs-$VERSION.tbz2` where `$VERSION` is the current version of OSPREY.
   This task is useful when building release versions of OSPREY.

 * **`archiveReleases`**:\
   This task uploads all the releases you currently have in `/build/releases`, including documentation releases,
   to the [DLab release archive][dlab-release-archive], if they're not there already.

 * **`downloadDocReleases`**:\
   This task downloads all the documentation releases from the [DLab release archive][dlab-release-archive].
   This task is useful for building the documentation website, which requires all the documentation releases.

 * **`buildWebsite`**:\
   This task builds the entire documentation website from the current contents of `/doc` and the
   documentation releases in `/build/releases`. This task will also generate the dynamic parts of the website, like
    * Download links to release files
    * Code documentation from the source file comments
   
   After building, the finished website files will appear in `/build/website-release`.

 * **`deployWebsite`**:\
   This task copies the documentation website in `/build/website-release` to the [DLab filesystem][dlab-docsite].

 * **`hugoServer`**:\
   This task starts a local Hugo development server, for quick feedback when editing Markdown documents for the
   website. Once started, the server stays running until you [tell it to stop][gradle-stop-task].

   {{% notice note %}}
   Dynamically generated content like the code documentation or download links will not be available
   on the Hugo server. Only static content from the `/doc` folder will be available.
   {{% /notice %}}

[dlab-release-archive]: {{< ref "/contributing/dlab-filesystem#release-archive" >}}
[dlab-docsite]: {{< ref "/contributing/dlab-filesystem#docsite" >}}


## Updating the live documentation website

To update the live documentation website, follow these steps:

 * Make your changes to the Markdown documents in `/doc/content`.
 * Review your changes locally by running the Hugo server with the `hugoServer` Gradle task.
 * Dowload any previous documentation releases with the `downloadDocReleases` Gradle task.
 * Build the documentation website with the `buildWebsite` Gradle task.
 * Upload the documentation website with the `deployWebsite` Gradle task.

Note that these tasks will reqiure access to the [DLab filesystem][dlab-filesystem].

[dlab-filesystem]: {{< ref "/contributing/dlab-filesystem" >}}
