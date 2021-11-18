
import sys
import re

import dataclasses
import typing as t

# https://niklasrosenstein.github.io/docspec/specification/
import docspec

from pydoc_markdown.interfaces import Processor, Resolver

import osprey_pydoc.javadoc


@dataclasses.dataclass
class OspreyProcessor(Processor):

	_tag = re.compile(r'\$\{([^\}]*)\}')

	def __init__(self):
		self._vtable = {
			'class_javadoc': self._class_javadoc,
			'type_java': self._type_java,
			'field_javadoc': self._field_javadoc,
			'arg_field_javadoc': self._arg_field_javadoc,
			'args_fields_javadoc': self._args_fields_javadoc,
			'returns_method_java': self._returns_method_java,
			'method_javadoc': self._method_javadoc,
			'arg_java': self._arg_java,
			'args_java': self._args_java,
			'arg_javadoc': self._arg_javadoc,
			'default': self._default
		}


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
			raise Exception('unrecognized macro function %s' % func_name)

		# parse the function arguments
		# allow nested list syntax, eg
		#   arg1
		#   arg1, arg2
		#   [arg1a, arg1b], arg2
		#   arg1(still_arg1,still_arg1), arg2
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
				next = []
				args_stack[-1].append(next)
				args_stack.append(next)
			elif c == ']':
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

		class_path = javadoc.Path(args[0])

		# lookup the class in the javadoc
		c = javadoc.get_class(class_path)
		if c is None:
			raise Exception('unknown java class: %s' % class_path)

		return c.javadoc

	def _type_java(self, args):

		# render the type with links to other docs where possible
		class_path = javadoc.Path(args[0])

		# lookup the class
		c = javadoc.get_class(class_path)
		if c is None:
			raise Exception('unknown java class: %s' % class_path)

		return _render_type(c.type)


	def _field_javadoc(self, args):

		field_path = javadoc.Path(args[0])

		# lookup the field in the javadoc
		field = javadoc.get_field_or_throw(field_path)

		return field.javadoc


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
		if arg.default_value == 'useJavaDefault':
			if field.initializer is not None:
				arg.default_value = _render_value(field.initializer)
			else:
				arg.default_value = 'None'

		# render the type
		if arg_type is not None:
			rendered_type = _literal_type(arg_type)
		else:
			rendered_type = _render_type(field.type)

		return '%s %s: %s' % (arg_name, rendered_type, field.javadoc)


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

		return _render_type(method.returns)


	def _method_javadoc(self, args):

		method = javadoc.get_method_or_throw(javadoc.Path(args[0]))

		if method.javadoc is None:
			raise Exception('java method %s has no javadoc' % method_path)

		return method.javadoc.header


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
			rendered_type = _render_type(arg.type)

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


	def _arg_javadoc(self, args):

		method = javadoc.get_method_or_throw(javadoc.Path(args[0]))
		arg_name = args[1]

		# lookup the arg javadoc
		try:
			arg = method.javadoc.params[arg_name]
		except (KeyError, AttributeError):
			raise Exception("can't find arg %s in java method %s\n\tavailable args: %s" % (arg_name, method_path, method.javadoc.params.keys()))

		return arg.description


	def _default(self, args):

		arg_name = args[0]
		value = args[1]

		# find the argument in the current node
		candidates = [a for a in self._current_node.args if a.name == arg_name]
		if len(candidates) <= 0:
			raise Exception("can't find arg %s in python funciton %s" % (arg_name, self._current_node.name))

		candidates[0].default_value = value


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
_URL_PREFIX = '../../java'

def _render_type(type):

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
		out.append('[%s](%s/%s)' % (type_name, _URL_PREFIX, type.url))
	else:
		out.append('`%s`' % type_name)

	# render the params
	if type.params is not None:
		out.append('<')
		for param in type.params:
			out.extend(_render_type(param))
		out.append('>')

	return ''.join(out)


def _literal_type(type):
	return '`%s`' % type


def _render_value(type):
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
