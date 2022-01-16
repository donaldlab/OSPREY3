
import sys
import re

import dataclasses
import typing as t

from nr.databind.json import JsonSerializer

# https://niklasrosenstein.github.io/docspec/specification/
import docspec

from pydoc_markdown.interfaces import Processor, Resolver

import osprey_pydoc.javadoc
import osprey_pydoc.kdoc


# pydoc_markdown apparently uses some ultra hacky library to emulate interfaces in python
# python doesn't provide interfaces natively, of course, because it's a freaking DUCK-TYPED LANGUAGE!!
# these nr.interfaces/databind libraries are apparently ultra-complicated, so OF COURSE they stopped working on our plugin code
# no idea why, but the best workaround I've found so far is to ignore all the garbage bloatware and
# Just Quack Like a Duck
# ie, just add process() and init() methods to our Processor class and be done with it
# but first, we have to make their configuration deserializer recognize our class by adding a custom deserializer
def _deserialize(mapper, node):
	# pydoc_markdown hides error information here for Some Unfathomable Reason
	# so explicitly capture and display it before letting pydoc_markdown waste our time
	try:
		return OspreyProcessor()
	except Exception as ex:
		print('exception!', ex, file=sys.stderr)
		raise Exception('OspreyProcessor instantiation failed')


@JsonSerializer(deserialize=_deserialize)
class OspreyProcessor:
#class OspreyProcessor(Processor):
	# NOTE: don't try to "implement" the Processor "interface"
	# here there be dragons!
	# the "interface" implementation is appa+rently very brittle and prone to breakage

	_tag = re.compile(r'\$\{([^\}]*)\}')

	def __init__(self):
		self._vtable = {
			'class_javadoc': self._class_javadoc,
			'class_kdoc': self._class_kdoc,
			'type_java': self._type_java,
			'type_jre': self._type_jre,
			'type_kotlin': self._type_kotlin,
			'type_kotlin_stdlib': self._type_kotlin_stdlib,
			'field_javadoc': self._field_javadoc,
			'arg_field_javadoc': self._arg_field_javadoc,
			'args_fields_javadoc': self._args_fields_javadoc,
			'returns_method_java': self._returns_method_java,
			'returns_func_kotlin': self._returns_func_kotlin,
			'method_javadoc': self._method_javadoc,
			'func_kdoc': self._func_kdoc,
			'arg_java': self._arg_java,
			'args_java': self._args_java,
			'arg_kotlin': self._arg_kotlin,
			'args_kotlin': self._args_kotlin,
			'receiver_kotlin': self._receiver_kotlin,
			'arg_javadoc': self._arg_javadoc,
			'enum_java': self._enum_java,
			'enum_kotlin': self._enum_kotlin,
			'default': self._default,
			'prop_kotlin': self._prop_kotlin,
			'link_func_kotlin': self._link_func_kotlin
		}

	def init(self, context):
		pass

	def process(self, modules: t.List[docspec.Module], resolver: t.Optional[Resolver]) -> None:
		docspec.visit(modules, self._process)


	def _process(self, node: docspec.ApiObject) -> bool:

		if node.docstring is None:
			return

		self._current_node = node

		# replace all expressions in the docstring, eg ${expr}
		node.docstring = self._tag.sub(self._expr, node.docstring)

		self._current_node = None


	def _expr(self, match):

		expr = match.group(1).strip()

		# parse the function name
		pivot = expr.find('(')
		if pivot <= 0 or expr[-1] != ')':
			raise Exception('expression is not a function call: %s' % expr)
		func_name = expr[0:pivot]

		# look up the function
		try:
			func = self._vtable[func_name]
		except KeyError:
			raise Exception('unrecognized macro function `%s`, try one of:\n\t%s' % (func_name, '\n\t'.join(self._vtable.keys())))

		# parse the function arguments
		# allow nested list syntax, eg
		#   arg1
		#   arg1, arg2
		#   [arg1a, arg1b], arg2
		#   arg1(still_arg1,still_arg1), arg2
		#   arg1(still_arg1[],still_arg1), arg2
		func_args = expr[pivot + 1:-1]
		args = []
		args_stack = [args]
		in_arg = False
		in_paren = False
		for c in func_args:
			if c == ',':
				if in_paren:
					args_stack[-1][-1] += c
				else:
					in_arg = False
			elif c == '(':
				in_paren = True
				args_stack[-1][-1] += c
			elif c == ')':
				in_paren = False
				args_stack[-1][-1] += c
			elif c == '[':
				if in_paren:
					args_stack[-1][-1] += c
				else:
					next = []
					args_stack[-1].append(next)
					args_stack.append(next)
			elif c == ']':
				if in_paren:
					args_stack[-1][-1] += c
				else:
					args_stack.pop()
					in_arg = False
			elif c == '\n' or c == '\t':
				pass # ignore some kinds of whitespace
			elif c == ' ':
				if in_arg:
					args_stack[-1][-1] += c
			else:
				if in_arg:
					args_stack[-1][-1] += c
				else:
					# start a new arg
					args_stack[-1].append(str(c))
					in_arg = True


		# actually call the function with the parsed arguments
		return func(args)

		# not a recognized expression
		raise Exception('unrecognized expression: %s' % match.group(0))


	def _class_javadoc(self, args):

		c = javadoc.get_class_or_throw(javadoc.Path(args[0]))

		return _render_javadoc(c.javadoc)


	def _class_kdoc(self, args):

		c = kdoc.get_class_or_throw(kdoc.Id(args[0]))

		return _render_kdoc(c.kdoc)


	def _type_java(self, args):

		c = javadoc.get_class_or_throw(javadoc.Path(args[0]))

		return _render_type_java(c.type)


	def _type_jre(self, args):

		type_name = args[0]

		type = javadoc.Type({
			'name': type_name,
			'url': None # TODO: link to Oracle docs online?
		})

		return _render_type_java(type)


	def _type_kotlin(self, args):

		c = kdoc.get_class_or_throw(kdoc.Id(args[0]))

		return _render_type_kotlin(c.type)


	def _type_kotlin_stdlib(self, args):

		id = kdoc.Id(args[0])

		type = kdoc.Type({
			'name': id.raw,
			'url': None # TODO: link to Kotlin docs online? The url transformations are non-trivial though =(
		})

		return _render_type_kotlin(type)


	def _field_javadoc(self, args):

		field_path = javadoc.Path(args[0])

		# lookup the field in the javadoc
		field = javadoc.get_field_or_throw(field_path)

		return _render_javadoc(field.javadoc)


	def _arg_field_javadoc(self, args):

		# output something like
		# arg (int): The first argument.

		(pos_args, named_args) = _split_args(args)
		arg_name = pos_args[0]
		field_path = javadoc.Path(pos_args[1])

		# parse the named args
		arg_type = None
		for (name, value) in named_args:
			if name == 'type':
				arg_type = value

		# find the corresponding argument in the function
		try:
			arg = [arg for arg in self._current_node.args if arg.name == arg_name][0]
		except IndexError:
			raise Exception('unknown argument %s in python function %s' % (arg_name, self._current_node.name))

		# lookup the field in the javadoc
		field = javadoc.get_field_or_throw(field_path)

		# use the initializer as an argument default value
		if arg.default_value in ['useJavaDefault', '_useJavaDefault']:
			if field.initializer is not None:
				arg.default_value = _render_value_java(field.initializer)
			else:
				arg.default_value = 'None'

		# render the type
		if arg_type is not None:
			rendered_type = _literal_type(arg_type)
		else:
			rendered_type = _render_type_java(field.type)

		return '%s %s: %s' % (arg_name, rendered_type, _render_javadoc(field.javadoc))


	def _args_fields_javadoc(self, args):

		classname = args[0]
		args = args[1:]

		out = []
		for field_args in args:

			(pos_args, named_args) = _split_args(field_args)

			# deconstruct the argument pair
			if len(pos_args) == 1:
				arg_name = pos_args[0]
				field_name = pos_args[0]
			elif len(pos_args) == 2:
				arg_name = pos_args[0]
				field_name = pos_args[1]
			else:
				raise Exception('invalid field arguments: %s' % field_args)

			field_path = '%s#%s' % (classname, field_name)
			out.append(self._arg_field_javadoc(_join_args([arg_name, field_path], named_args)))

		return '\n'.join(out)


	def _returns_method_java(self, args):

		path = javadoc.Path(args[0])
		method = javadoc.get_method_or_throw(path)

		if method.returns is None:
			raise Exception("can't write return type, method %s has no return type (maybe it's a constructor?)" % path)

		return _render_type_java(method.returns)


	def _returns_func_kotlin(self, args):

		id = kdoc.Id(args[0])
		func = kdoc.get_func_or_throw(id)

		return _render_type_kotlin(func.returns)


	def _method_javadoc(self, args):

		method_path = javadoc.Path(args[0])
		method = javadoc.get_method_or_throw(method_path)

		if method.javadoc is None:
			raise Exception('java method %s has no javadoc' % method_path)

		return _render_javadoc(method.javadoc)


	def _func_kdoc(self, args):

		func_id = kdoc.Id(args[0])
		func = kdoc.get_func_or_throw(func_id)

		if func.kdoc is None:
			raise Exception('kotlin function %s has no kdoc' % func_id)

		return _render_kdoc(func.kdoc)


	def _arg_java(self, args):

		(pos_args, named_args) = _split_args(args)
		py_arg_name = pos_args[0]
		method = javadoc.get_method_or_throw(javadoc.Path(pos_args[1]))
		java_arg_name = pos_args[0]
		if len(args) >= 3:
			java_arg_name = pos_args[2]

		# parse the named args
		arg_type = None
		for (name, value) in named_args:
			if name == 'type':
				arg_type = value

		# lookup the arg in the java method
		arg = method.find_arg_or_throw(java_arg_name)

		# render the type
		if arg_type is not None:
			rendered_type = _literal_type(arg_type)
		else:
			rendered_type = _render_type_java(arg.type)

		# lookup the arg javadoc, if any
		try:
			arg_javadoc = method.javadoc.params[java_arg_name]
		except (KeyError, AttributeError):
			arg_javadoc = ''

		return '%s %s: %s' % (py_arg_name, rendered_type, arg_javadoc)


	def _args_java(self, args):

		path = javadoc.Path(args[0])
		method = javadoc.get_method_or_throw(path)
		args = args[1:]

		out = []
		for arg_args in args:
			out.append(self._arg_java([arg_args[0]] + [path.path] + arg_args[1:]))

		return '\n'.join(out)


	def _arg_kotlin(self, args):

		(pos_args, named_args) = _split_args(args)
		py_arg_name = pos_args[0]
		func = kdoc.get_func_or_throw(kdoc.Id(pos_args[1]))
		kotlin_arg_name = pos_args[0]
		if len(args) >= 3:
			kotlin_arg_name = pos_args[2]

		# parse the named args
		arg_type = None
		for (name, value) in named_args:
			if name == 'type':
				arg_type = value

		# lookup the arg in the kotlin function
		arg = func.find_arg_or_throw(kotlin_arg_name)

		# render the type
		if arg_type is not None:
			rendered_type = _literal_type(arg_type)
		else:
			rendered_type = _render_type_kotlin(arg.type)

		# lookup the arg kdoc, if any
		try:
			arg_kdoc = func.kdoc.params[kotlin_arg_name]
		except (KeyError, AttributeError):
			arg_kdoc = ''

		return '%s %s: %s' % (py_arg_name, rendered_type, arg_kdoc)


	def _args_kotlin(self, args):

		id = kdoc.Id(args[0])
		func = kdoc.get_func_or_throw(id)
		args = args[1:]

		out = []
		for arg_args in args:
			out.append(self._arg_kotlin([arg_args[0]] + [id.raw] + arg_args[1:]))

		return '\n'.join(out)


	def _receiver_kotlin(self, args):

		py_arg_name = args[0]
		func = kdoc.get_func_or_throw(kdoc.Id(args[1]))

		return '%s %s' % (py_arg_name, _render_type_kotlin(func.receiver))


	def _arg_javadoc(self, args):

		method = javadoc.get_method_or_throw(javadoc.Path(args[0]))
		arg_name = args[1]

		# lookup the arg javadoc
		try:
			arg = method.javadoc.params[arg_name]
		except (KeyError, AttributeError):
			raise Exception("can't find arg %s in java method %s\n\tavailable args: %s" % (arg_name, method_path, method.javadoc.params.keys()))

		return arg.description


	def _enum_java(self, args):

		path = javadoc.Path(args[0])
		c = javadoc.get_class_or_throw(path)

		out = ''

		# show the javadoc, if any
		if c.javadoc is not None:
			out += _render_javadoc(c.javadoc)
			out += '\n\n'

		# also show the enum values
		out += '**Enumeration:** %s' % _render_type_java(c.type)
		for field_name in c.field_names:

			# filter down to just the enum values
			field = c.field_or_throw(field_name)
			if field.type.name != c.type.name:
				continue

			# render the enum values as a markdown list
			out += '\n\n * **%s**:' % field_name
			if field is not None and field.javadoc is not None:
				out += ' '
				out += _render_javadoc(field.javadoc)

		return out


	def _enum_kotlin(self, args):

		id = kdoc.Id(args[0])
		c = kdoc.get_class_or_throw(id)

		out = ''

		# show the kdoc, if any
		if c.kdoc is not None:
			out += _render_kdoc(c.kdoc)
			out += '\n\n'

		# also show the enum values
		out += '**Enumeration:** %s' % _render_type_kotlin(c.type)
		for entry_name in c.entry_names:
			entry = c.entry_or_throw(entry_name)

			# render the enum values as a markdown list
			out += '\n\n * **%s**:' % entry.name
			if entry.kdoc is not None:
				out += ' '
				out += _render_kdoc(entry.kdoc)

		return out


	def _default(self, args):

		arg_name = args[0]
		value = args[1]

		# find the argument in the current node
		candidates = [a for a in self._current_node.args if a.name == arg_name]
		if len(candidates) <= 0:
			raise Exception("can't find arg %s in python funciton %s" % (arg_name, self._current_node.name))

		candidates[0].default_value = value


	def _prop_kotlin(self, args):

		id = kdoc.Id(args[0])
		prop = kdoc.get_prop_or_throw(id)

		return _render_type_kotlin(prop.type)


	def _link_func_kotlin(self, args):

		id = kdoc.Id(args[0])
		func = kdoc.get_func_or_throw(id)

		return '[%s](%s/%s)' % (func.name, _URL_PREFIX_KDOC, func.url)


def _split_args(args):

	# pull out the named optional arguments from the positional arguments
	named_args = [a.split('=') for a in args if '=' in a]
	pos_args = [a for a in args if '=' not in a]

	# check the named args
	for parts in named_args:
		if len(parts) != 2:
			raise Exception('invalid named argument: %s' % parts)

	return (pos_args, named_args)


def _join_args(pos_args, named_args):
	return pos_args + ['='.join(a) for a in named_args]


# use relative URLs here, not absoulte URLs, so the docs folders are copyable
_URL_PREFIX_JAVADOC = '../../java'
_URL_PREFIX_KDOC = '../../kotlin'

def _render_type_java(type):

	if type is None:
		raise Exception("can't render type info, no type")

	out = []

	# perform any java->python type transformations that are automatically applied by JPype
	try:
		type_name = {
			'java.lang.String': 'str',
			'boolean': 'bool',
			'double': 'float',
			'long': 'int'
		}[type.name]
	except KeyError:
		type_name = type.name

	# render the type name, and a link if possible
	if type.url is not None:
		out.append('[%s](%s/%s)' % (type_name, _URL_PREFIX_JAVADOC, type.url))
	else:
		out.append('`%s`' % type_name)

	# render the params
	if type.params is not None:
		out.append('<')
		for param in type.params:
			out.extend(_render_type_java(param))
		out.append('>')

	return ''.join(out)


def _render_type_kotlin(type):

	if type is None:
		raise Exception("can't render type info, no type")

	out = []

	# perform any java->python type transformations that are automatically applied by JPype
	try:
		type_name = {
			'java.lang.String': 'str',
			'kotlin.String': 'str',
			'kotlin.Boolean': 'bool',
			'kotlin.Float': 'float',
			'kotlin.Double': 'float',
			'kotlin.Int': 'int',
			'kotlin.Long': 'int'
		}[type.name]
	except KeyError:
		type_name = type.name

	# render the variance, if any
	if type.variance is not None:
		out.append('%s ' % type.variance)

	# render the type name, and a link if possible
	if type.url is not None:
		out.append('[%s](%s/%s)' % (type_name, _URL_PREFIX_KDOC, type.url))
	else:
		out.append('`%s`' % type_name)

	# render the params
	if type.params is not None:
		out.append('<')
		for param in type.params:
			out.extend(_render_type_kotlin(param))
		out.append('>')

	if type.nullable:
		out.append('?')

	return ''.join(out)


def _literal_type(type):
	return '`%s`' % type


def _render_value_java(type):
	# translate from Java to Python
	try:
		return {
			'null': 'None',
			'false': 'False',
			'true': 'True',
			'Double.POSITIVE_INFINITY': "float('inf')",
			'Double.NEGATIVE_INFINITY': "float('-inf')"
		}[type]
	except KeyError:
		return type


def _render_javadoc(javadoc):

	if javadoc is None:
		return ''

	# start with the raw javadoc text
	text = javadoc.text

	# replace citations
	for citation in javadoc.citations:
		text = text.replace(citation.text, _render_citation(citation.lines), 1)

	# replace warnings
	for warning in javadoc.warnings:
		text = text.replace(warning.text, _render_warning(warning.content), 1)

	# replace notes
	for note in javadoc.notes:
		text = text.replace(note.text, _render_note(note.content), 1)

	# replace todos
	for todo in javadoc.todos:
		text = text.replace(todo.text, _render_todo(todo.content), 1)

	# replace links
	for link in javadoc.links:
		if link.signature is not None:
			url = '%s#%s' % (link.type.url, link.signature)
		else:
			url = link.type.url
		repl = '[%s](%s/%s)' % (link.label, _URL_PREFIX_JAVADOC, url)
		text = text.replace(link.text, repl, 1)

	# add block elements to the end
	for warning in javadoc.warnings:
		if warning.is_block:
			text += '\n\n'
			text += _render_warning(warning.content)
	for note in javadoc.notes:
		if note.is_block:
			text += '\n\n'
			text += _render_note(note.content)
	for todo in javadoc.todos:
		if todo.is_block:
			text += '\n\n'
			text += _render_todo(todo.content)

	return text


def _render_kdoc(kdoc):

	if kdoc is None:
		return ''

	# TODO: do we need anything fancy here?

	return kdoc.text


def _render_notice(kind, content):
	return r'{{% notice ' + kind + r' %}}' + '\n' + content + '\n' + r'{{% /notice %}}'


def _render_citation(lines):
	# NOTE: add a \ between lines to preserve line breaks
	return _render_notice('info', '\\\n'.join(lines))


def _render_warning(content):
	return _render_notice('warning', content)


def _render_note(content):
	return _render_notice('note', content)


def _render_todo(content):
	return _render_notice('tip', '**TODO** ' + content)
