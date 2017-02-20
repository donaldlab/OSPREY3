
import os
import docutils, sphinx
import javalang
import re


ast_cache = {}
doctags = []

# java->python type translations, see:
# http://jpype.readthedocs.io/en/latest/userguide.html#type-conversion
java_to_python_types = {
	'byte': 'int',
	'short': 'int',
	'int': 'int',
	'long': 'long',
	'float': 'float',
	'double': 'float',
	'char': 'str',
	'String': 'str',
	'boolean': 'bool'
}

java_to_python_constants = {
	'true': 'True',
	'false': 'False',
	'null': 'None',
	'Double.POSITIVE_INFINITY': "float('inf')",
	'Double.NEGATIVE_INFINITY': "-float('inf')"
}


# Sphinx application API
# www.sphinx-doc.org/en/1.5.1/extdev/appapi.html

def setup(app):

	# configs
	app.add_config_value('javadoc_sources_dir', '.', 'env')
	app.add_config_value('javadoc_package_prefix', '', 'env')

	app.add_domain(JavadocDomain)

	app.add_role('java:classdoc', ClassdocRole())
	app.add_role('java:fielddoc', FielddocRole())
	app.add_role('java:methoddoc', MethoddocRole())

	app.connect('autodoc-process-signature', autodoc_signature_handler)
	app.connect('autodoc-process-docstring', autodoc_docstring_handler)

	# builder doctags should come before default!
	# (because builder doctags emit default doctags)
	doctags.append(BuilderOptionDoctag())
	doctags.append(BuilderReturnDoctag())
	doctags.append(DefaultDoctag())

	return { 'version': '0.1' }


# make sure warnings have the correct source info in autodoc events
class AutodocHandlerWarn():
	
	def __init__(self, app, func):
		self.app = app
		self.func = func

	
	def __call__(self, msg, lineno=None, cause=None):
		modname = self.func.__module__
		srcpath = sphinx.pycode.ModuleAnalyzer.for_module(modname).srcname
		loc = '%s:docstring of %s.%s' % (srcpath, modname, self.func.__name__)
		if lineno is not None:
			loc = '%s:%d' % (loc, lineno)
		if cause is not None:
			msg += '\nCause: ' + str(cause)
		self.app.warn(msg, location=loc)
		
		

def autodoc_signature_handler(app, what, name, obj, options, signature, return_annotation):

	# get the function arg info
	try:
		argspec = sphinx.util.inspect.getargspec(obj)
	except TypeError:
		# not a function, don't mess with the signature
		return

	# convert argspec to a mutable format (argspecs are not editable by default)
	argspec = argspec._asdict()
	if argspec['defaults'] is not None:
		argspec['defaults'] = list(argspec['defaults'])

	warn = AutodocHandlerWarn(app, obj)

	# run the docstring handlers to expand doctags before modifying the signature
	doclines = obj.__doc__.split('\n')
	autodoc_docstring_handler(app, what, name, obj, options, doclines, for_signature=True)

	# process each doctag in order
	for tag in doctags:

		# look for tags in the docstring, like:
		# :default argname: default value
		for i in range(len(doclines)):
			line = doclines[i]

			try:
				name, args, text, indent, sublines = Doctag.parse_lines(doclines, i)
				if name != tag.name:
					raise ValueError
			except (KeyError, ValueError):
				# not a doctag, just keep going
				continue

			# curry the warn function for this line
			def curried_warn(msg, cause=None):
				warn(msg, lineno=i, cause=cause)

			tag.handle_signature(app, args, text, obj, argspec, curried_warn)
				
	# build the new signature from the modified argspec
	return (sphinx.ext.autodoc.formatargspec(obj, *argspec.values()), None)


def autodoc_docstring_handler(app, what, name, obj, options, doclines, for_signature=False):

	warn = AutodocHandlerWarn(app, obj)

	# process each doctag in order
	for tag in doctags:

		newlines = []
		skip_lines = 0

		# look for tags in the docstring, like:
		# :default argname: default value
		for i in range(len(doclines)):

			# skip lines if needed
			if skip_lines > 0:
				skip_lines -= 1
				continue

			line = doclines[i]

			try:
				name, args, text, indent, sublines = Doctag.parse_lines(doclines, i)
				if name != tag.name:
					raise ValueError
			except (KeyError, ValueError):
				# not our doctag, just keep going
				newlines.append(line)
				continue

			# curry the warn function for this line
			def curried_warn(msg, cause=None):
				warn(msg, lineno=i, cause=cause)

			# what does the tag want to do with this line?
			result = tag.handle_docstring(app, args, text, indent, line, sublines, curried_warn, for_signature)
			if result is not None:

				# add the result from the doctag
				if isinstance(result, str):
					newlines.append(result)
				else:
					newlines.extend(result)

				# skip the sublines
				skip_lines = len(sublines)


		# replace the existing lines
		doclines[:] = newlines

	# DEBUG: for docstring transformations
	#if for_signature is False and obj.__name__ == 'TemplateLibrary':
	#	print('DOCSTRING')
	#	for line in doclines:
	#		print(line)


def resolve_typename(typename, package, imports, config):

	# is this a primitive type?
	try:
		return java_to_python_types[typename]
	except KeyError:
		# not a primitive
		pass
	
	# try a java class
	ref = ImportResolver(package, imports, config).resolve(typename)
	return ':java:ref:`%s`' % ref.full_classname


class Doctag():

	@staticmethod
	def parse_lines(lines, lineno):

		# analyze the first line for doctag params
		line = lines[lineno]

		def get_indent(line):
			try:
				return re.findall(r'^\s+', line)[0]
			except IndexError:
				return ''

		indent = get_indent(line)
		line = line.strip()

		if not line.startswith(':'):
			raise ValueError

		line_parts = line.split(':')
		if len(line_parts) < 3:
			raise ValueError

		name_and_args = line_parts[1].strip().split(' ')
		name = name_and_args[0]
		args = name_and_args[1:]
		text = ':'.join(line_parts[2:]).strip()

		# get sublines (have greater indent)
		last_sublineno = lineno
		for i in range(lineno + 1, len(lines)):
			nextline = lines[i]
			last_sublineno = i
			if len(nextline.strip()) == 0:
				continue
			# did we hit a line with same or lesser indent?
			# NOTE: this won't work with mixed tabs and spaces,
			# but let's just hope that never happens
			if len(get_indent(nextline)) <= len(indent):
				last_sublineno = i
				break

		sublines = lines[lineno + 1:last_sublineno]

		# if all sublines are empty, remove them
		nonempty_sublines = [l for l in sublines if len(l.strip()) > 0]
		if len(nonempty_sublines) == 0:
			sublines = []

		return name, args, text, indent, sublines


	def handle_signature(self, app, args, text, func, argspec, warn):
		# override me if you want
		pass

	
	def handle_docstring(self, app, args, text, indent, line, sublines, warn, for_signature=False):
		# do nothing to remove this line from the docstring
		# or override me if you want
		pass


class DefaultDoctag(Doctag):

	name = 'default'

	def map_constants(self, val):
		
		# apply common java->python transformations for well-known constants
		if val in java_to_python_constants:
			val = java_to_python_constants[val]

		return val
	

	def resolve(self, target, config):

		# get the java field
		ref = JavaRef(expand_classname(target, config))
		ast = get_class_ast(ref, config)
		field = ast.find_field(ref.membername)

		# get the default value from the field declaration/initialization
		init = field.declarator.initializer
		if isinstance(init, javalang.tree.Literal):

			# initialized to literal value, use directly
			return self.map_constants(init.value)

		elif isinstance(init, javalang.tree.MemberReference):

			# initialized to some expression like foo or foo.bar
			if init.qualifier is '':
				val = init.member
			else:
				val = '%s.%s' % (init.qualifier, init.member)

			return self.map_constants(val)

		else:
			raise ValueError("field initialized with a %s, don't know what to do" % init)


	def set_default(self, argspec, name, val):

		# write defaults without quoting strings
		class Literal():
			def __init__(self, val):
				self.val = val
			def __repr__(self):
				return self.val

		names = argspec['args']
		defaults = argspec['defaults']

		offset = len(names) - len(defaults)

		for i in range(len(defaults)):
			if names[offset + i] == name:
				defaults[i] = Literal(val)
				return

		raise KeyError
	

	def handle_signature(self, app, args, text, func, argspec, warn):

		defaults = argspec['defaults']
		if defaults is None:
			return

		# resolve any references in the first line of the text
		ref_regex = re.compile(r':java:default:`([^`]+)`')
		def ref_resolver(match):
			try:
				target = match.group(1)
				return self.resolve(target, app.config)
			except (KeyError, ValueError, FileNotFoundError) as e:
				warn("can't resolve java ref: %s" % target, cause=e)
				return match.group(0)
		text = ref_regex.sub(ref_resolver, text)

		# apply default value to the signature
		try:
			self.set_default(argspec, args[0], text)
		except KeyError:
			warn("can't find argument '%s'" % args[0])
			return


	def handle_docstring(self, app, args, text, indent, line, sublines, warn, for_signature=False):

		# for signatures, leave our line in, so the signature doctags can find them
		if for_signature:
			return [line] + sublines

		# otherwise, remove them
		else:
			return None


class BuilderOptionDoctag(Doctag):

	name = 'builder_option'

	def handle_docstring(self, app, args, text, indent, line, sublines, warn, for_signature=False):

		if len(args) < 2:
			warn('Not enough args for builder option. Need python arg name, builder field ref')
			return None

		# parse the args
		argname = args[0]
		ref = JavaRef(expand_classname(args[1], app.config))

		# look up the builder class
		try:
			ast = get_class_ast(ref, app.config)
		except (KeyError, ValueError, FileNotFoundError) as e:
			warn("can't find builder class: %s" % ref.full_classname, cause=e)
			return

		# lookup the builder field
		try:
			field = ast.find_field(ref.membername)
		except KeyError:
			warn("can't find field %s" % ref)
			return

		rst = []

		if for_signature:

			# if the field has an assignment, grab the value for the arg default
			if field.declarator.initializer is not None:
				rst.append('%s:default %s: :java:default:`%s`' % (indent, argname, ref))

		else:

			# use the javadoc as the param desc, if available
			if field.documentation is None:
				javadoc_desc = ''
			else:
				javadoc = Javadoc(field.documentation, ast.package, ast.imports, app.config)
				javadoc_desc = javadoc.description
			rst.append('%s:param %s: %s' % (indent, argname, javadoc_desc))

			# add sublines after the :param: tag
			rst.extend(sublines)

			# use the field type as the argument type
			try:
				typename = field.type.name
				typename = resolve_typename(typename, ast.package, ast.imports, app.config)
				rst.append('%s:type %s: %s' % (indent, argname, typename))
			except ValueError as e:
				warn("can't resolve java type '%s' referenced in class '%s'" % (typename, ref), cause=e)
				return

		return rst


class BuilderReturnDoctag(Doctag):

	name = 'builder_return'

	def handle_docstring(self, app, args, text, indent, line, sublines, warn, for_signature=False):

		# look up the builder class
		ref = JavaRef(expand_classname(args[0], app.config))
		try:
			ast = get_class_ast(ref, app.config)
		except (KeyError, ValueError, FileNotFoundError) as e:
			warn("can't find builder class: %s" % ref, cause=e)
			return

		# get the build method
		try:
			buildmethod = ast.find_method('build')
			typename = buildmethod.return_type.name
			typename = resolve_typename(typename, ast.package, ast.imports, app.config)
			return '%s:rtype: %s' % (indent, typename)
		except KeyError:
			warn("can't find build method on builder class: %s" % ref)
		except ValueError as e:
			warn("can't resolve java type '%s' referenced in class '%s'" % (typename, ref), cause=e)


def expand_classname(classname, config):
	
	# expand the classname if we need
	if classname[0] == '.':
		classname = config.javadoc_package_prefix + classname

	return classname


class JavaRef():
	
	def __init__(self, target, membername=None):

		# accepts javadoc style targets, eg:
		# package.Class
		# package.Class#member
		# package.Class$Inner
		# package.Class$Inner$Inner#member

		# 'full classnames' are fully qualified, eg:
		# package.Class
		# package.Class$Inner

		# 'classnames' have no packages, but have nesting chains, eg:
		# Class
		# Class$Inner$Inner

		# 'simple classnames' are just the name of the class, eg:
		# Class
		# Inner

		# split off the package
		parts = target.split('.')
		if len(parts) > 1:
			self.package = '.'.join(parts[:-1])
			target = parts[-1]
		else:
			self.package = None

		# get the member name
		if membername is not None:
			self.membername = membername
		else:
			parts = target.split('#')
			if len(parts) > 1:
				self.membername = parts[1]
				target = parts[0]
			else:
				self.membername = None

		self.classname = target

		# parse the nested classes
		self.simple_classnames = self.classname.split('$')
		self.outer_classname = self.simple_classnames[0]
		self.simple_classname = self.simple_classnames[-1]
		self.simple_inner_classnames = self.simple_classnames[1:]

		# get the fully-qualified classnames
		if self.package is None:
			self.full_classname = None
			self.full_outer_classname = None
		else:
			self.full_classname = '.'.join([self.package, self.classname])
			self.full_outer_classname = '.'.join([self.package, self.outer_classname])
	

	def __str__(self):
		out = self.classname
		if self.package is not None:
			out = '%s.%s' % (self.package, out)
		if self.membername is not None:
			out = '%s#%s' % (out, self.membername)
		return out


	def unprefixed_full_classname(self, prefix):
		if self.full_classname.startswith(prefix):
			return self.full_classname[len(prefix) + 1:]
		else:
			return self.full_classname

	
	def source_file(self):
		return '%s.java' % self.full_outer_classname.replace('.', '/')


class ImportResolver():

	def __init__(self, package, imports, config):
		self.package = package
		self.imports = imports
		self.config = config


	def has_source(self, ref):
		path = os.path.join(self.config.javadoc_sources_dir, ref.source_file())
		return os.path.exists(path)

	
	def resolve(self, target):

		ref = JavaRef(target)

		# find an import that resolves the target
		for imp in self.imports:
			if imp.path.endswith('.%s' % ref.simple_classname):
				return JavaRef(imp.path, membername=ref.membername)

		# no matching import, check package
		ref = JavaRef('%s.%s' % (self.package.name, ref.classname), membername=ref.membername)
		if self.has_source(ref):
			return ref
		else:
			raise ValueError("can't resolve java class against package or imports: %s" % target)


def read_file(path):
	with open(path, 'r') as file:
		return file.read()


def get_class_ast(ref, config):

	# convert ref to a filepath
	path = os.path.join(config.javadoc_sources_dir, ref.source_file())

	# read syntax tree from the cache, or from the file
	try:
		ast = ast_cache[path]
	except KeyError:
		ast = javalang.parse.parse(read_file(path))
		ast_cache[path] = ast

	# save package and imports for later
	package = ast.package
	imports = ast.imports

	def find_type(name, classtypes):
		for classtype in classtypes:
			if classtype.name == name:
				return classtype
		raise KeyError
	
	# get the root type in the compilation unit
	try:
		ast = find_type(ref.outer_classname, ast.types)
	except KeyError:
		raise KeyError("can't find outer class %s in source file %s" % (ref.outer_classname, ref.source_file()))
	
	# get nested types if needed
	for name in ref.simple_inner_classnames:
		subtypes = [member for member in ast.body if isinstance(member, javalang.tree.TypeDeclaration)]
		try:
			ast = find_type(name, subtypes)
		except KeyError:
			raise KeyError("can't find inner class %s in outer class %s in source file %s" % (name, ast.name, ref.source_file()))
	
	# reference the package and imports in the returned type
	ast.package = package
	ast.imports = imports

	# add some helper methods
	def find_field(name):
		for field in ast.fields:
			for declarator in field.declarators:
				if declarator.name == name:

					# attach the matched declarator to the field
					field.declarator = declarator

					return field
		raise KeyError("can't find field: %s" % name)
	ast.find_field = find_field

	def find_method(name):
		found_methods = []
		for method in ast.methods:
			if method.name == name:
				found_methods.append(method)
		if len(found_methods) == 0:
			raise KeyError("can't find method: %s" % name)
		elif len(found_methods) > 1:
			raise KeyError("multiple method overloads found for: %s, can't resolve" % name)
		else:
			return found_methods[0]
	ast.find_method = find_method

	return ast


class Javadoc():

	def __init__(self, text, package, imports, config):

		if text is None:
			raise ValueError('javadoc cannot be None')

		self.text = text
		self.import_resolver = ImportResolver(package, imports, config)
		self.parsed = javalang.javadoc.parse(text)


	@property
	def description(self):
		return self._translate(self.parsed.description)


	def _translate(self, text):

		# translate links from javadoc to sphinx
		link_regex = re.compile(r'\{@link ([^\}]+)\}')
		def link_resolver(match):
			
			# what's being linked to?
			target = match.group(1)

			# resolve against imports to get a full ref
			ref = self.import_resolver.resolve(target)

			# build the rst
			rst = []
			rst.append(':java:ref:`%s`' % str(ref))
			return '\n'.join(rst)

		return link_regex.sub(link_resolver, text)


def is_public(thing):
	return 'public' in thing.modifiers


def should_show(thing):
	return thing.documentation is not None and is_public(thing)


def make_method_signature(method):
	return method.name


def parse_rst(rst, settings):
	
	# if rst is a list of strings, join it
	if not isinstance(rst, str):
		rst = '\n'.join(rst)

	# parse the rst and return the nodes
	docSource = 'I come from the water'
	doc = docutils.utils.new_document(docSource, settings)
	docutils.parsers.rst.Parser().parse(rst, doc)

	# if there's just a single paragraph node at the root of the doc, unwrap it
	# so our parsed rst can be placed inline with existing rst
	if len(doc.children) == 1 and isinstance(doc.children[0], docutils.nodes.paragraph):
		return doc.children[0].children
	
	return doc.children


class ParsingDirective(docutils.parsers.rst.Directive):

	@property
	def env(self):
		return self.state.document.settings.env


	@property
	def config(self):
		return self.env.config


	def parse(self, rst):
		return parse_rst(rst, self.state.document.settings)


class JavaClassDirective(ParsingDirective):
	
	has_content = True

	def run(self):

		ref = JavaRef(expand_classname(self.content[0], self.config))
		ast = get_class_ast(ref, self.config)

		rst = []

		# class name header
		rst.append(ast.name)
		rst.append('=' * len(ast.name))

		rst.append('')

		showedSomething = False

		# show fields
		fields = [field for field in ast.fields if should_show(field)]
		if len(fields) > 0:

			# show fields
			rst.append('Properties')
			rst.append('----------')
			rst.append('')
			showedSomething = True

			for field in fields:

				# show the field name and javadoc
				for decl in field.declarators:
					rst.append('.. py:attribute:: %s' % decl.name)
					rst.append('')
					if field.documentation is not None:
						javadoc = Javadoc(field.documentation, ast.package, ast.imports, self.config)
						rst.append('\t' + javadoc.description)
						rst.append('')

		# show methods
		methods = [method for method in ast.methods if should_show(method)]
		if len(methods) > 0:

			rst.append('Methods')
			rst.append('-------')
			rst.append('')
			showedSomething = True

			for method in methods:

				# show the method signature and javadoc
				rst.append('.. py:method:: %s' % make_method_signature(method))
				rst.append('')
				if method.documentation is not None:
					javadoc = Javadoc(method.documentation, ast.package, ast.imports, self.config)
					rst.append('\t' + javadoc.description)
					rst.append('')

		# show enum constants
		if isinstance(ast, javalang.tree.EnumDeclaration):

			rst.append('Constants')
			rst.append('-----------')
			rst.append('')
			showedSomething = True

			for value in ast.body.constants:
				
				# show the value name and javadoc
				rst.append('.. py:attribute:: %s' % value.name)
				rst.append('')
				if value.documentation is not None:
					javadoc = Javadoc(value.documentation, ast.package, ast.imports, self.config)
					rst.append('\t' + javadoc.description)
					rst.append('')

		if not showedSomething:

			rst.append('*(This topic does not yet have documentation)*')

		return self.parse(rst)


# how to write function roles:
# http://docutils.sourceforge.net/docs/howto/rst-roles.html

class JavaRole():

	def __call__(self, name, rawtext, text, lineno, inliner, options={}, content=[]):

		settings = inliner.document.settings

		# add some convenience attributes for subclasses
		self.config = settings.env.config
	
		# add some convenience methods for subclasses
		def warn(msg, cause=None):

			# format the warning for rst
			rst = '*(%s)*' % msg

			# add cause info for the console
			if cause is not None:
				msg += '\n\tcause: ' + str(cause)

			inliner.reporter.warning(msg, line=lineno)

			return rst

		self.warn = warn

		# parse the rst and return docutils nodes
		return parse_rst(self.make_rst(text), settings), []


	def make_rst(self, text):
		raise Exception('implement me for %s' % self.__class__.__name__)


class ClassdocRole(JavaRole):

	def make_rst(self, text):

		# find the class
		try:
			ref = JavaRef(expand_classname(text, self.config))
			ast = get_class_ast(ref, self.config)
		except (KeyError, ValueError, FileNotFoundError) as e:
			return self.warn("can't find java class: %s" % text, cause=e)
		
		# look for the javadoc
		if ast.documentation is None:
			return self.warn("class %s has no javadoc" % ref)

		# parse the javadoc
		try:
			return Javadoc(ast.documentation, ast.package, ast.imports, self.config).description
		except (KeyError, ValueError) as e:
			return self.warn("can't parse javadoc for class: %s" % ref, cause=e)


class FielddocRole(JavaRole):

	def make_rst(self, text):

		# find the field
		try:
			ref = JavaRef(expand_classname(text, self.config))
			ast = get_class_ast(ref, self.config)
			field = ast.find_field(ref.membername)
		except (KeyError, ValueError, FileNotFoundError) as e:
			return self.warn("can't find field: %s" % text, cause=e)
		
		# look for the javadoc
		if field.documentation is None:
			return self.warn("field %s has no javadoc" % ref)

		# parse the javadoc
		try:
			return Javadoc(field.documentation, ast.package, ast.imports, self.config).description
		except (KeyError, ValueError) as e:
			return self.warn("can't parse javadoc for field: %s" % ref, cause=e)


class MethoddocRole(JavaRole):

	def make_rst(self, text):

		# find the method
		try:
			ref = JavaRef(expand_classname(text, self.config))
			ast = get_class_ast(ref, self.config)
			method = ast.find_method(ref.membername)
		except (KeyError, ValueError, FileNotFoundError) as e:
			return self.warn("can't find method: %s" % text, cause=e)

		# look for the javadoc
		if method.documentation is None:
			return self.warn("method %s has no javadoc" % ref)

		# parse the javadoc
		try:
			return Javadoc(method.documentation, ast.package, ast.imports, self.config).description
		except (KeyError, ValueError) as e:
			return self.warn("can't parse javadoc for method: %s" % ref, cause=e)


class RefRole(sphinx.roles.XRefRole):

	def __call__(self, typ, rawtext, text, lineno, inliner, options={}, content=[]):

		# get the reporter from the inliner
		# so the warnings have the correct source info when using autodoc
		self.reporter = inliner.reporter

		# then make a convenience method for warnings
		def warn(msg):
			self.reporter.warning(msg, line=lineno)
		self.warn = warn

		return sphinx.roles.XRefRole.__call__(self, typ, rawtext, text, lineno, inliner, options, content)


	def process_link(self, env, refnode, has_explicit_title, title, target):

		# resolve the target
		ref = JavaRef(expand_classname(target, env.config))
		resolved = ResolvedXref()
		resolved.docpath = 'api.' + ref.unprefixed_full_classname(env.config.javadoc_package_prefix)
		resolved.docpath = resolved.docpath.replace('$', '.')

		if ref.membername is not None:
			resolved.text = ref.membername
			resolved.anchor = ref.membername
		else:
			resolved.text = ref.simple_classname

		# see if the docpath exists
		path = env.doc2path(resolved.docpath)
		if not os.path.exists(path):
			self.warn('cross-reference does not exist: %s' % path)

		# attach the resolution to the node so we can find it again in Domain.resolve_xref()
		refnode.resolved = resolved

		return title, target
	

class ResolvedXref():

	def __init__(self):
		self.docpath = None
		self.text = None
		self.anchor = ''
		self.title = None


# Sphinx domain API:
# http://www.sphinx-doc.org/en/1.5.1/extdev/domainapi.html

class JavadocDomain(sphinx.domains.Domain):

	name = 'java'
	label = 'Osprey javadoc processor'

	directives = {
		'class': JavaClassDirective
	}

	roles = {
		'ref': RefRole()
	}


	def resolve_xref(self, env, fromdocname, builder, reftype, target, node, contnode):

		# get the ref resolution from our JavaRole
		resolved = node.resolved

		# set the link text in the node thingy
		contnode.replace(
			contnode.children[0],
			docutils.nodes.Text(resolved.text)
		)

		# return the ref node
		return sphinx.util.nodes.make_refnode(
			builder,
			fromdocname,
			resolved.docpath,
			resolved.anchor,
			contnode,
			resolved.title
		)
 
