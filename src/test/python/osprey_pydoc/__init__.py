
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
			'arg_field_javadoc': self._arg_field_javadoc,
			'args_fields_javadoc': self._args_fields_javadoc,
			'returns_method_java': self._returns_method_java,
			'method_javadoc': self._method_javadoc,
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


	def _arg_field_javadoc(self, args):

		# output something like
		# arg (int): The first argument.

		arg_name = args[0]
		field_path = javadoc.Path(args[1])

		# find the corresponding argument in the function
		try:
			arg = [arg for arg in self._current_node.args if arg.name == arg_name][0]
		except IndexError:
			raise Exception('unknown argument %s in python function %s' % (arg_name, self._current_node.name))

		# lookup the field in the javadoc
		field = javadoc.get_field(field_path)
		if field is None:
			raise Exception('unknown java field: %s' % field_path)

		# does the field have an initializer?
		if field.initializer is not None:

			# use the initializer as an argument default value
			if arg.default_value == 'UseJavaDefault':
				arg.default_value = field.initializer

		return '%s %s: %s' % (arg_name, _render_type(field.type), field.javadoc)


	def _args_fields_javadoc(self, args):

		classname = args[0]
		args = args[1:]

		out = []
		for pair in args:

			# deconstruct the argument pair
			if len(pair) == 1:
				arg_name = pair[0]
				field_name = pair[0]
			elif len(pair) == 2:
				arg_name = pair[0]
				field_name = pair[1]
			else:
				raise Exception('invaild argument description: %s' % pair)

			field_path = '%s#%s' % (classname, field_name)
			out.append(self._arg_field_javadoc([arg_name, field_path]))

		return '\n'.join(out)


	def _returns_method_java(self, args):

		method = javadoc.get_method_or_throw(javadoc.Path(args[0]))

		return _render_type(method.returns)


	def _method_javadoc(self, args):

		method = javadoc.get_method_or_throw(javadoc.Path(args[0]))

		if method.javadoc is None:
			raise Exception('java method %s has no javadoc' % method_path)

		return method.javadoc.header


	def _arg_javadoc(self, args):

		method = javadoc.get_method_or_throw(javadoc.Path(args[0]))
		arg_name = args[1]

		# lookup the arg javadoc
		try:
			arg = method.javadoc.params[arg_name]
		except KeyError:
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


# use relative URLs here, not absoulte URLs, so the docs folders are copyable
_URL_PREFIX = '../../java'

def _render_type(type):

	out = []

	# render the type name, and a link if possible
	if type.url is not None:
		out.append('[%s](%s/%s)' % (type.name, _URL_PREFIX, type.url))
	else:
		out.append('`%s`' % type.name)

	# render the params
	if type.params is not None:
		out.append('<')
		for param in type.params:
			out.extend(_render_type(param))
		out.append('>')

	return ''.join(out)
