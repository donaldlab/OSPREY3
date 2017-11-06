
import os
import docutils, sphinx
import javalang
import re
import six


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

	# for module variables, look for special DocstringXXX instances
	if what == 'data':
		if obj.__class__.__name__ == 'DocstringJavaDefault':
			options['annotation'] = ' = %s' % DefaultDoctag().resolve(repr(obj), app.config)
			return

	# otherwise, get the function arg info
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
				if isinstance(result, six.string_types):
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


def resolve_type(typeast, rootast, config):

	# java wildcards in generics can cause a null type here
	if typeast is None:
		return "?"

	name = typeast.name

	# is this a primitive type?
	try:
		name = java_to_python_types[typeast.name]
	except KeyError:
		# don't have a mapping, keep the raw type name
		pass
	
	# is it a class?
	if isinstance(typeast, javalang.tree.ReferenceType):
		try:
			ref = ImportResolver(rootast, config).resolve(typeast.name)
			name = ':java:ref:`%s`' % ref.full_classname
		except ValueError:
			# can't resolve name just keep raw type
			name = '``%s``' % name

		# handle type params
		if typeast.arguments is not None:
			type_args = [resolve_type(arg.type, rootast, config) for arg in typeast.arguments]
			name = '%s < %s >' % (name, ' , '.join(type_args))

	# handle arrays
	for i in range(len(typeast.dimensions)):
		# NOTE: need a space before the [] to prevent RST parser errors
		name = name + ' []'

	return name


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


	def render_ast_expression(self, ast):
		
		if isinstance(ast, javalang.tree.Literal):

			# expr is literal value
			return ast.value

		elif isinstance(ast, javalang.tree.MemberReference):

			# expr is foo or foo.bar
			val = ast.member
			if len(ast.qualifier) > 0:
				val = '%s.%s' % (ast.qualifier, val)
			return val

		elif isinstance(ast, javalang.tree.ClassCreator):
			
			# expr is new Foo()
			return '%s()' % ast.type.name

		elif isinstance(ast, javalang.tree.MethodInvocation):

			# expr is foo() or foo.bar() or foo.bar(arg1, arg2)
			args = []
			for arg in ast.arguments:
				args.append(self.render_ast_expression(arg))
			val = '%s(%s)' % (ast.member, ', '.join(args))
			if len(ast.qualifier) > 0:
				val = '%s.%s' % (ast.qualifier, val)
			return val

		elif isinstance(ast, javalang.tree.BinaryOperation):

			# expr is e.g., 4 + 4

			left = self.render_ast_expression(ast.operandl)
			if isinstance(ast.operandl, javalang.tree.BinaryOperation):
				left = '(%s)' % left

			right = self.render_ast_expression(ast.operandr)
			if isinstance(ast.operandr, javalang.tree.BinaryOperation):
				right = '(%s)' % right

			return '%s %s %s' % (left, ast.operator, right)

		else:
			raise ValueError("java expression is a %s, don't know how to render" % ast)
	

	def resolve(self, target, config):

		# get the java field
		ref = JavaRef(expand_classname(target, config))
		ast = get_class_ast(ref, config)
		field = ast.find_field(ref.membername)

		# render the initializer expression
		return self.map_constants(self.render_ast_expression(field.declarator.initializer))


	def set_default(self, argspec, name, val):

		# write defaults without quoting strings
		class Literal():
			def __init__(self, val):
				self.val = val
			def __repr__(self):
				return self.val

		names = argspec['args']
		defaults = argspec['defaults']

		num_args = len(names) - len(defaults)
		num_kwargs = len(defaults)

		# is this an arg?
		for i in range(num_args):
			if names[i] == name:
				# do nothing, can't set defaults for regular args
				return

		# is this a keyword arg?
		for i in range(num_kwargs):
			if names[num_args + i] == name:
				# get a kwarg, set the default
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
		except javalang.parser.JavaSyntaxError as e:
			warn("can't parse java source: %s" % ref, cause=e)
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
				desc = ''
			else:
				desc = ':java:fielddoc:`%s`' % ref
			rst.append('%s:param %s: %s' % (indent, argname, desc))

			# add sublines after the :param: tag
			rst.extend(sublines)

			# use the field type as the argument type
			typename = resolve_type(field.type, ast, app.config)
			rst.append('%s:type %s: %s' % (indent, argname, typename))

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
		except javalang.parser.JavaSyntaxError as e:
			warn("can't parse java source: %s" % ref, cause=e)
			return

		# get the build method
		try:
			buildmethod = ast.find_method('build')
			typename = resolve_type(buildmethod.return_type, ast, app.config)
			return '%s:rtype: %s' % (indent, typename)
		except KeyError:
			warn("can't find build method on builder class: %s" % ref)


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

		# method arguments and generic types are not handled though

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

		# remove any method args from the member name
		if self.membername is not None:
			self.membername = self.membername.split('(')[0]

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

	def __init__(self, ast, config):
		self.ast = ast
		self.config = config


	def has_source(self, ref):
		path = os.path.join(self.config.javadoc_sources_dir, ref.source_file())
		return os.path.exists(path)

	
	def find_inner_type(self, ast, simple_name):

		try:
			return ast.find_inner_class(simple_name)
		except KeyError:
			# nope, not there, check recursively
			for inner_type in ast.inner_types:
				try:
					return self.find_inner_type(inner_type, simple_name)
				except KeyError:
					continue

		# didn't find anything
		raise KeyError

	
	def resolve(self, target):

		ref = JavaRef(target)

		# first, check the main type in the ast
		if self.ast.name == ref.simple_classname:
			return JavaRef(self.ast.ref.full_classname, membername=ref.membername)

		# then, check ast inner-types
		try:
			innerref = self.find_inner_type(self.ast, ref.simple_classname).ref
			return JavaRef(innerref.full_classname, membername=ref.membername)
		except KeyError:
			# can't find it
			pass

		# then, check ast outer-types
		try:
			t = self.ast.outer_type
			try:
				outerref = self.find_inner_type(t, ref.simple_classname).ref
				return JavaRef(outerref.full_classname, membername=ref.membername)
			except KeyError:
				# not there, try another outer
				t = t.outer_type
				pass

		except AttributeError:
			# no more outer types
			pass

		# find an import that resolves the target
		for imp in self.ast.imports:
			if imp.path.endswith('.%s' % ref.simple_classname):

				# java import statements don't separate inner classes from packages,
				# so try to detect them here
				check_ref = JavaRef(imp.path, membername=ref.membername)
				while not self.has_source(check_ref):

					# is there a . to replace with a $?
					if '.' not in check_ref.package:
						# nope, can't resolve this import
						break

					# assume one more class is an inner class and try again
					new_classname = '%s$%s' % (check_ref.package, check_ref.classname)
					check_ref = JavaRef(new_classname, membername=ref.membername)

				return check_ref

		# no matching import, check package
		ref = JavaRef('%s.%s' % (self.ast.package.name, ref.classname), membername=ref.membername)
		if self.has_source(ref):
			return ref
		else:
			raise ValueError("can't resolve java class against package or imports: %s" % target)


def read_file(path):
	with open(path, 'r') as file:
		return file.read()


# augment javalang types with helper methods

def _ast_find_type(name, classtypes):
	for classtype in classtypes:
		if classtype.name == name:
			return classtype
	raise KeyError


@property
def _ast_inner_types(self):
	return [member for member in self.body if isinstance(member, javalang.tree.TypeDeclaration)]
javalang.tree.TypeDeclaration.inner_types = _ast_inner_types


def _ast_find_inner_class(self, name):
	inner_class = _ast_find_type(name, self.inner_types)

	# copy over metadata
	inner_class.package = self.package
	inner_class.imports = self.imports
	inner_class.ref = JavaRef('%s$%s' % (self.ref.full_classname, name))

	return inner_class
javalang.tree.TypeDeclaration.find_inner_class = _ast_find_inner_class


def _ast_find_method(self, name, require_javadoc=False):
	found_methods = []
	for method in self.methods:
		if require_javadoc and method.documentation is None:
			continue
		if method.name == name:
			found_methods.append(method)
	if len(found_methods) == 0:
		raise KeyError("can't find method: %s" % name)
	elif len(found_methods) > 1:
		raise KeyError("multiple method overloads found for: %s, can't resolve" % name)
	else:
		return found_methods[0]
javalang.tree.TypeDeclaration.find_method = _ast_find_method


def _ast_find_field(self, name):
	for field in self.fields:
		for declarator in field.declarators:
			if declarator.name == name:

				# attach the matched declarator to the field
				field.declarator = declarator

				return field
	raise KeyError("can't find field: %s" % name)
javalang.tree.TypeDeclaration.find_field = _ast_find_field


def get_class_ast(ref, config):

	# convert ref to a filepath
	path = os.path.join(config.javadoc_sources_dir, ref.source_file())

	# read syntax tree from the cache, or from the file
	try:
		ast = ast_cache[path]
	except KeyError:
		ast = javalang.parse.parse(read_file(path))
		ast_cache[path] = ast

	# get the root type in the compilation unit
	try:
		typeast = _ast_find_type(ref.outer_classname, ast.types)
	except KeyError:
		raise KeyError("can't find outer class %s in source file %s" % (ref.outer_classname, ref.source_file()))

	# copy package and imports to the type
	typeast.package = ast.package
	typeast.imports = ast.imports

	# add a ref to the type
	typeast.ref = JavaRef(ref.full_outer_classname)

	# get nested types if needed
	for name in ref.simple_inner_classnames:
		try:
			inner_typeast = typeast.find_inner_class(name)
			inner_typeast.outer_type = typeast
			typeast = inner_typeast
		except KeyError:
			raise KeyError("can't find inner class %s in outer class %s in source file %s" % (name, typeast.name, ref.source_file()))
	
	# add some helper methods to ast classes
	return typeast


class Javadoc():

	def __init__(self, text, ast, config):

		if text is None:
			raise ValueError('javadoc cannot be None')

		self.text = text
		self.import_resolver = ImportResolver(ast, config)
		self.parsed = javalang.javadoc.parse(text)

		# start with the main javadoc text
		desc = [self._translate(self.parsed.description)]

		# handle see tags
		try:
			for text in self.parsed.tags['see']:
				desc.append('See %s' % self._resolve_link(text))

		except KeyError:
			pass

		# translate any notes, warnings, etc into rst
		for tagname in ['note', 'warning', 'todo']:
			try:
				for text in self.parsed.tags[tagname]:
					text = text.replace('\n', ' ')
					desc.append('')
					desc.append('.. %s:: %s' % (tagname, self._translate(text)))
			except KeyError:
				pass

		self.description = '\n'.join(desc)

	
	def arg(self, name):
		for argname, argdesc in self.parsed.params:
			if argname == name:
				return argdesc
		raise KeyError


	def _resolve_link(self, target):

		# add enclosing class info if needed
		if target.startswith('#'):
			ref = JavaRef(self.import_resolver.ast.ref.full_outer_classname, membername=target[1:])
		
		# otherwise, resolve against imports to get a full ref
		else:
			ref = self.import_resolver.resolve(target)

		# build the rst
		return ':java:ref:`%s`' % str(ref)
		
	
	def _resolve_all(self, text, tagname, resolver):
		regex = re.compile(r'\{@%s ([^\}]+)\}' % tagname)
		def resolver_helper(match):
			return resolver(match.group(1))
		return regex.sub(resolver_helper, text)
		

	def _translate(self, text):

		# escape * characters
		text = text.replace(r'*', r'\*')

		# translate links from javadoc to sphinx
		text = self._resolve_all(text, 'link', self._resolve_link)

		# translate citations to rst
		citations = {}
		def cite_resolver(key_and_citation):

			# parse the citation
			parts = key_and_citation.split(' ')
			key = parts[0]
			citation = ' '.join(parts[1:]).replace('\n', ' ')

			# save the citaiton for later, so we can append it
			citations[key] = citation

			# put just the key inline
			return '[%s]_' % key

		text = self._resolve_all(text, 'cite', cite_resolver)

		# append the citations to the end
		for key, citation in citations.items():
			text += '\n\n.. [%s] %s' % (key, citation)
	
		return text


def is_public(thing, ast):
	if isinstance(ast, javalang.tree.ClassDeclaration):
		return 'public' in thing.modifiers
	elif isinstance(ast, javalang.tree.InterfaceDeclaration):
		return 'private' not in thing.modifiers and 'protected' not in thing.modifiers


def should_show(thing, ast):
	return is_public(thing, ast) # and thing.documentation is not None


def is_constant(thing, ast):
	return 'static' in thing.modifiers and 'final' in thing.modifiers


def parse_rst(rst, rst_source, settings):
	
	# if rst is a list of strings, join it
	if not isinstance(rst, six.string_types):
		rst = '\n'.join(rst)

	# DEBUG: for investigating warnings in dynamically-generated RST
	#print('DYNAMIC RST: %s' % rst_source)
	#lineno = 1
	#for line in rst.split('\n'):
	#	print('%3d  %s' % (lineno, line))
	#	lineno += 1

	# parse the rst and return the nodes
	docSource = 'dynamically-generated-rst-' + rst_source
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


	def parse(self, rst, rst_source):
		return parse_rst(rst, rst_source, self.state.document.settings)


	def warn(self, msg, cause=None):
		if cause is not None:
			msg = '%s\nCause: %s' % (msg, str(cause))
		self.state.memo.reporter.warning(msg)
		return self.parse(msg, 'warn')


class JavaClassDirective(ParsingDirective):
	
	has_content = True

	def run(self):

		try:
			ref = JavaRef(expand_classname(self.content[0], self.config))
			ast = get_class_ast(ref, self.config)
		except (KeyError, ValueError, IOError) as e:
			return self.warn("can't find class: %s" % ref, cause=e)
		except javalang.parser.JavaSyntaxError as e:
			return self.warn("can't parse java source: %s" % ref, cause=e)

		rst = []

		# class name header
		rst.append(ast.name)
		rst.append('=' * len(ast.name))
		rst.append('')

		# show the full classname
		rst.append('Java class: ``%s``' % ref)
		rst.append('')

		# show the class javadoc, if any
		if ast.documentation is not None:
			try:
				javadoc = Javadoc(ast.documentation, ast, self.config)
				if javadoc is not None:
					rst.append(self.indent(self.format(javadoc.description)))
					rst.append('')
			except ValueError as e:
				self.warn("Can't parse javadoc for class %s" % ref, cause=e)

		# add nested content
		if len(self.content) > 1:
			rst.extend(self.content[1:])
			rst.append('')

		rstlen = len(rst)

		# show constructors
		constructors = [constructor for constructor in ast.constructors if should_show(constructor, ast)]
		if len(constructors) > 0:
			self.show_header(rst, 'Constructors')
			
			# add a link to the tutorial for how to call java constructors from python
			rst.append('*This is a Java class which has been exposed to the Python API. For more information')
			rst.append('on how to call constructors on Java objects, see the tutorial:* :ref:`constructors`.')
			rst.append('')
			rst.append('-'*10)
			rst.append('')

			for constructor in constructors:
				self.show_method(rst, ref, constructor, ast)

		# show constants
		constants = [field for field in ast.fields if should_show(field, ast) and is_constant(field, ast)]
		if len(constants) > 0:
			self.show_header(rst, 'Constants')
			for constant in constants:
				self.show_field_constant(rst, ref, constant, ast)

		# show fields
		fields = [field for field in ast.fields if should_show(field, ast) and not is_constant(field, ast)]
		if len(fields) > 0:
			self.show_header(rst, 'Properties')
			for field in fields:
				self.show_field(rst, ref, field, ast)

		# show methods
		methods = [method for method in ast.methods if should_show(method, ast)]
		if len(methods) > 0:
			self.show_header(rst, 'Methods')
			for method in methods:
				self.show_method(rst, ref, method, ast)

		# show enum constants
		if isinstance(ast, javalang.tree.EnumDeclaration):
			self.show_header(rst, 'Constants')
			for value in ast.body.constants:
				self.show_enum_constant(rst, ref, value, ast)
				

		# add a message if we didn't show anything
		if len(rst) == rstlen:
			rst.append('*(This topic does not yet have documentation)*')
			rst.append('')

		# DEBUG
		#if ref.classname.startswith('SimpleConfSpace'):
		#	print('CLASSDOC')
		#	for line in rst:
		#		print(line)

		return self.parse(rst, 'class-%s' % ref)
	
		
	def show_header(self, rst, title):
		rst.append(title)
		rst.append('_' * len(title))
		rst.append('')


	def indent(self, rst, num=1):

		prefix = '\t'*num

		if isinstance(rst, six.string_types):
			return prefix + rst.replace('\n', '\n' + prefix)
		else:
			return [prefix + line for line in rst]


	def format(self, javadoc):

		# parse preformatted blocks
		pre_regex = re.compile(r'<pre>(.+?)</pre>', re.DOTALL | re.IGNORECASE)
		def parse_pre(match):
			text = match.group(1)
			return '::\n' + self.indent(text)
		javadoc = pre_regex.sub(parse_pre, javadoc)

		return javadoc
	
	
	def show_field(self, rst, ref, field, ast):

		# get the javadoc, if any
		javadoc = None
		if field.documentation is not None:
			try:
				javadoc = Javadoc(field.documentation, ast, self.config)
			except ValueError as e:
				for decl in field.declarators:
					self.warn("Can't parse javadoc for field %s#%s" % (ref, decl.name), cause=e)

		# show the field name and javadoc
		for decl in field.declarators:
			rst.append('.. py:attribute:: %s.%s' % (ref.classname, decl.name))
			rst.append('')
			if javadoc is not None:
				rst.append(self.indent(self.format(javadoc.description)))
				rst.append('')

			rst.append('\t:type: %s' % resolve_type(field.type, ast, self.config))
			rst.append('')


	def show_field_constant(self, rst, ref, field, ast):

		# get the javadoc, if any
		javadoc = None
		if field.documentation is not None:
			try:
				javadoc = Javadoc(field.documentation, ast, self.config)
			except ValueError as e:
				self.warn("Can't parse javadoc for field %s#%s" % (ref, field.name), cause=e)

		# show the field name and javadoc
		for decl in field.declarators:

			# render the value
			renderer = DefaultDoctag()
			value = renderer.map_constants(renderer.render_ast_expression(decl.initializer))

			rst.append('.. py:attribute:: %s.%s = %s' % (ref.classname, decl.name, value))
			rst.append('')
			if javadoc is not None:
				rst.append(self.indent(self.format(javadoc.description)))
				rst.append('')

			rst.append('\t:type: %s' % resolve_type(field.type, ast, self.config))
			rst.append('')


	def show_method(self, rst, ref, method, ast):

		# get the javadoc, if any
		javadoc = None
		if method.documentation is not None:
			try:
				javadoc = Javadoc(method.documentation, ast, self.config)
			except ValueError as e:
				self.warn("Can't parse javadoc for method %s#%s" % (ref, method.name), cause=e)

		# show the method signature
		args = [arg.name for arg in method.parameters]
		sig = '%s.%s(%s)' % (ref.classname, method.name, ', '.join(args))
		rst.append('.. py:method:: %s' % sig)
		rst.append('')

		# show the javadoc desc, if any
		if javadoc is not None:
			rst.append(self.indent(self.format(javadoc.description)))
			rst.append('')

		# method arg descriptions, types
		for arg in method.parameters:
			try:
				desc = self.format(javadoc.arg(arg.name))
				rst.append('\t:param %s: %s' % (arg.name, desc))
			except (AttributeError, KeyError):
				# no javadoc for this arg
				rst.append('\t:param %s:' % arg.name)

			rst.append('\t:type %s: %s' % (arg.name, resolve_type(arg.type, ast, self.config)))

		# return type
		if isinstance(method, javalang.tree.MethodDeclaration) and method.return_type is not None:
			rst.append('\t:rtype: %s' % resolve_type(method.return_type, ast, self.config))

		rst.append('')


	def show_enum_constant(self, rst, ref, value, ast):
		
		# get the javadoc, if any
		javadoc = None
		if value.documentation is not None:
			try:
				javadoc = Javadoc(value.documentation, ast, self.config)
			except ValueError as e:
				self.warn("Can't parse javadoc for enum value %s.%s" % (ref, value.name), cause=e)

		# show the value name and javadoc
		rst.append('.. py:attribute:: %s.%s' % (ref.classname, value.name))
		rst.append('')
		if javadoc is not None:
			rst.append(self.indent(self.format(javadoc.description)))
			rst.append('')


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
		return parse_rst(self.make_rst(text), 'JavaRole-' + name, settings), []


	def make_rst(self, text):
		raise Exception('implement me for %s' % self.__class__.__name__)


class ClassdocRole(JavaRole):

	def make_rst(self, text):

		# find the class
		try:
			ref = JavaRef(expand_classname(text, self.config))
			ast = get_class_ast(ref, self.config)
		except (KeyError, ValueError, IOError) as e:
			return self.warn("can't find java class: %s" % text, cause=e)
		except javalang.parser.JavaSyntaxError as e:
			return self.warn("can't parse java source: %s" % ref, cause=e)
		
		# look for the javadoc
		if ast.documentation is None:
			return self.warn("class %s has no javadoc" % ref)

		# parse the javadoc
		try:
			return Javadoc(ast.documentation, ast, self.config).description
		except (KeyError, ValueError) as e:
			return self.warn("can't parse javadoc for class: %s" % ref, cause=e)


class FielddocRole(JavaRole):

	def make_rst(self, text):

		# find the field
		try:
			ref = JavaRef(expand_classname(text, self.config))
			ast = get_class_ast(ref, self.config)
			field = ast.find_field(ref.membername)
		except (KeyError, ValueError, IOError) as e:
			return self.warn("can't find field: %s" % text, cause=e)
		except javalang.parser.JavaSyntaxError as e:
			return self.warn("can't parse java source: %s" % ref, cause=e)
		
		# look for the javadoc
		if field.documentation is None:
			return self.warn("field %s has no javadoc" % ref)

		# parse the javadoc
		try:
			return Javadoc(field.documentation, ast, self.config).description
		except (KeyError, ValueError) as e:
			return self.warn("can't parse javadoc for field: %s" % ref, cause=e)


class MethoddocRole(JavaRole):

	def make_rst(self, text):

		# find the method
		try:
			ref = JavaRef(expand_classname(text, self.config))
			ast = get_class_ast(ref, self.config)
			method = ast.find_method(ref.membername, require_javadoc=True)
		except (KeyError, ValueError, IOError) as e:
			return self.warn("can't find method: %s" % text, cause=e)
		except javalang.parser.JavaSyntaxError as e:
			return self.warn("can't parse java source: %s" % ref, cause=e)

		# look for the javadoc
		if method.documentation is None:
			return self.warn("method %s has no javadoc" % ref)

		# parse the javadoc
		try:
			return Javadoc(method.documentation, ast, self.config).description
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

		# is this one of our classes?
		if ref.package is not None and ref.package.startswith(env.config.javadoc_package_prefix):

			# resolve the link
			name = ref.unprefixed_full_classname(env.config.javadoc_package_prefix)
			filename = ('api.' + name).replace('$', '.')
			resolved.docpath = filename

			# see if the docpath exists, and create it if not
			filepath = env.doc2path(filename)
			if not os.path.exists(filepath):

				f = open(filepath, 'w')
				f.write("\n:orphan:\n\n.. java:class:: .%s\n\n" % name)
				f.close()

				self.warn('created RST file for class %s at %s\nrun `make html` again to generate HTML for this RST.' % (name, filepath))

		if ref.membername is not None:
			resolved.text = ref.membername
			resolved.anchor = '%s.%s' % (ref.simple_classname, ref.membername)
		else:
			resolved.text = ref.simple_classname

		resolved.title = str(ref)

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

		# if no resolution, just make a literal node
		if resolved.docpath is None:
			node = docutils.nodes.literal('', '', internal=True)
			node.append(contnode)
			return node

		# return a ref node
		return sphinx.util.nodes.make_refnode(
			builder,
			fromdocname,
			resolved.docpath,
			resolved.anchor,
			contnode,
			resolved.title
		)
 

# HACKHACK: add a convenience method to JavaSyntaxErrors so we can render them correctly
def java_syntax_error_str(self):
	return '%s at %s' % (self.description, self.at)
javalang.parser.JavaSyntaxError.__str__ = java_syntax_error_str


# HACKHACK: work around a missing feature in javalang:
# parsing of default interface methods

# see: https://github.com/c2nes/javalang/blob/master/javalang/parser.py

def javalang_parse_interface_method_declarator_rest(self):

	parameters = self.parse_formal_parameters()
	array_dimension = self.parse_array_dimension()
	throws = None

	if self.try_accept('throws'):
		throws = self.parse_qualified_identifier_list()

	body = None
	if self.would_accept('{'):
		body = self.parse_block()
	else:
		self.accept(';')

	return javalang.tree.MethodDeclaration(
		parameters=parameters,
		throws=throws,
		body=body,
		return_type=javalang.tree.Type(dimensions=array_dimension)
	)
javalang.parser.Parser.parse_interface_method_declarator_rest = javalang_parse_interface_method_declarator_rest

def javalang_parse_void_interface_method_declarator_rest(self):

	parameters = self.parse_formal_parameters()
	throws = None

	if self.try_accept('throws'):
		throws = self.parse_qualified_identifier_list()

	body = None
	if self.would_accept('{'):
		body = self.parse_block()
	else:
		self.accept(';')

	return javalang.tree.MethodDeclaration(
		parameters=parameters,
		throws=throws,
		body=body
	)
javalang.parser.Parser.parse_void_interface_method_declarator_rest = javalang_parse_void_interface_method_declarator_rest

