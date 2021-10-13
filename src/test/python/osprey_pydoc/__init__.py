
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

	tag = re.compile(r'\$\{([^\}]*)\}')
	tag_func = re.compile(r'(\w+)\(([^\)]*)\)')


	def process(self, modules: t.List[docspec.Module], resolver: t.Optional[Resolver]) -> None:
		docspec.visit(modules, self._process)


	def _process(self, node: docspec.ApiObject) -> bool:

		if node.docstring is None:
			return

		self._current_node = node

		# replace all expressions in the docstring, eg ${expr}
		node.docstring = self.tag.sub(self._expr, node.docstring)

		self._current_node = None


	def _expr(self, match):

		# try to match a function expression
		func = self.tag_func.match(match.group(1))
		if func is not None:
			func_name = func.group(1)

			# parse the function arguments
			# allow nested list syntax, eg
			#   arg1
			#   arg1, arg2
			#   [arg1a, arg1b], arg2
			func_args = func.group(2)
			args = []
			args_stack = [args]
			in_arg = False
			for c in func_args:
				if c == ',':
					in_arg = False
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
						args_stack[-1].append(str(c))
						in_arg = True

			return self._func(func_name, args)

		# not a recognized expression
		raise Exception('unrecognized expression: %s' % match.group(0))


	def _func(self, name, args):
		if name == 'class_javadoc':
			return self._class_javadoc(args)
		elif name == 'arg_field_javadoc':
			return self._arg_field_javadoc(args)
		elif name == 'args_fields_javadoc':
			return self._args_fields_javadoc(args)
		elif name == 'returns_java':
			return self._returns_java(args)
		elif name == 'returns_method_javadoc':
			return self._returns_method_javadoc(args)
		else:
			raise Exception('unrecognized macro function: %s' % name)


	def _class_javadoc(self, args):

		class_path = javadoc.Path(args[0])

		# lookup the class in the javadoc
		c = javadoc.get_class(class_path)
		if c is None:
			raise Exception('unknown java class: %s' % class_path)

		return c.javadoc


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

		return '%s (%s): %s' % (arg_name, field.type, field.javadoc)


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


	def _returns_method_javadoc(self, args):

		method_path = javadoc.Path(args[0])

		# lookup the method in the javadoc
		method = javadoc.get_method(method_path)
		if method is None:
			raise Exception('unknown java method: %s' % method_path)

		return '(%s):' % method.returns


	def _returns_java(self, args):

		class_path = javadoc.Path(args[0])

		# lookup the class in javadoc
		c = javadoc.get_class(class_path)
		if c is None:
			raise Exception('unknown java class: %s' % class_path)

		return '(%s):' % c.name