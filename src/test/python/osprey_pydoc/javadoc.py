
import pathlib
import json


# load the parsed javadoc
with open('build/doc/javadoc.json', 'r') as f:
	_javadoc = json.load(f)


default_package = 'edu.duke.cs.osprey'


class Path:

	def __init__(self, path: str):

		# path is something like, eg:
		#  package.Class
		#  p1.p2.Class
		#  p1.p2.Class$InnerClass
		#  package.Class#member
		#  package.Class$Inner1$Inner2#member
		#  .relative.package.Class

		# apply the default package root if needed
		if path[0] == '.':
			path = '%s.%s' % (default_package, path[1:])

		self.path = path

		# get the fully-qualified class name
		if '#' in path:
			parts = path.split('#')
			self.classname = parts[0]
			self.member = parts[-1]
		else:
			self.classname = path
			self.member = None


	def __str__(self):
		return self.path

	def __repr__(self):
		return 'Path[classname=%s, member=%s]' % (self.classname, self.member)


class Class:

	def __init__(self, json):
		self.type = Type(json['type'])
		try:
			self.javadoc = Javadoc(json['javadoc'])
		except KeyError:
			self.javadoc = None


def get_class(path: Path):
	try:
		c = _javadoc[path.classname]
	except KeyError:
		return None
	return Class(c)


def get_class_or_throw(path: Path):
	c = get_class(path)
	if c is None:
		raise Exception('unknown java class: %s' % path.classname)
	return c


class Field:

	def __init__(self, name, json):
		self.name = name
		self.type = Type(json['type'])
		try:
			self.javadoc = Javadoc(json['javadoc'])
		except KeyError:
			self.javadoc = None
		try:
			self.initializer = json['initializer']
		except KeyError:
			self.initializer = None


def get_field(path: Path):
	try:
		field = _javadoc[path.classname]['fields'][path.member]
	except KeyError:
		return None
	return Field(path.member, field)


def get_field_or_throw(path: Path):
	field = get_field(path)
	if field is None:
		raise Exception('unknown java field: %s' % path)
	return field


class Method:

	def __init__(self, id, json):
		self.id = id
		self.name = json['name']
		self.signature = json['signature']
		try:
			self.returns = Type(json['returns'])
		except KeyError:
			self.returns = None
		try:
			self.javadoc = Javadoc(json['javadoc'])
		except KeyError:
			self.javadoc = None
		self.args = [MethodArg(a) for a in json['args']]

	def find_arg(self, name):
		for arg in self.args:
			if arg.name == name:
				return arg
		return None

	def find_arg_or_throw(self, name):
		arg = self.find_arg(name)
		if arg is None:
			raise Exception('no argument named %s found in method %s\n\ttry one of: %s' % (name, self.id, [a.name for a in self.args]))
		return arg

class MethodArg:

	def __init__(self, json):
		self.name = json['name']
		self.type = Type(json['type'])


def get_methods(path: Path):

	try:
		c = _javadoc[path.classname]
	except KeyError:
		return []

	methods = c['methods']

	try:
		return [Method(path.member, methods[path.member])]
	except KeyError:

		# didn't get a direct hit, try to match just the method name without the type signature
		ids = [id for (id, m) in methods.items() if m['name'] == path.member]

		return [Method(id, methods[id]) for id in ids]


def get_method(path: Path):

	methods = get_methods(path)

	if len(methods) <= 0:
		return None
	elif len(methods) > 1:
		ids = [m.id for m in methods]
		raise Exception('multiple overloads found for method %s in %s\ntry one of:\n\t%s' % (path.member, path.classname, '\n\t'.join(ids)))

	return methods[0]


def get_method_or_throw(path: Path):

	# look for the class or throw, since get_method() won't throw on a missing class
	get_class_or_throw(path)

	method = get_method(path)
	if method is None:
		raise Exception('unknown java method: %s\ntry one of:\n\t%s' % (path, '\n\t'.join(get_method_ids(path))))
	return method


def get_method_ids(path: Path):

	try:
		c = _javadoc[path.classname]
	except KeyError:
		return []

	return c['methods'].keys()


class Type:

	def __init__(self, json):

		try:
			self.name = json['name']
		except (TypeError, KeyError):
			# no object here, the json bit must be the string name
			self.name = json

		try:
			self.url = json['url']
		except (KeyError, TypeError):
			self.url = None

		try:
			self.params = [Type(param) for param in json['params']]
		except (KeyError, TypeError):
			self.params = None


class Javadoc:

	def __init__(self, json):

		self.text = json['text']

		# read the params, if any
		try:
			params_json = json['params']
		except KeyError:
			params_json = []
		self.params = {}
		for param_json in params_json:
			param = Param(param_json)
			self.params[param.name] = param

		# read the citations, if any
		try:
			citations_json = json['citations']
		except KeyError:
			citations_json = []
		self.citations = []
		for citation_json in citations_json:
			self.citations.append(Citation(citation_json))

		# read the warnings, if any
		try:
			warnings_json = json['warnings']
		except KeyError:
			warnings_json = []
		self.warnings = []
		for warning_json in warnings_json:
			self.warnings.append(Warning(warning_json))

		# read the notes, if any
		try:
			notes_json = json['notes']
		except KeyError:
			notes_json = []
		self.notes = []
		for note_json in notes_json:
			self.notes.append(Note(note_json))

		# read the todos, if any
		try:
			todos_json = json['todos']
		except KeyError:
			todos_json = []
		self.todos = []
		for todo_json in todos_json:
			self.todos.append(Todo(todo_json))

		# read the links, if any
		try:
			links_json = json['links']
		except KeyError:
			links_json = []
		self.links = []
		for link_json in links_json:
			self.links.append(Link(link_json))


class Param:

	def __init__(self, json):
		self.text = json['text']
		self.name = json['name']
		self.description = json['description']


class Citation:

	def __init__(self, json):
		self.text = json['text']
		self.lines = json['lines']


class Warning:

	def __init__(self, json):
		self.text = json['text']
		self.content = json['content']


class Note:

	def __init__(self, json):
		self.text = json['text']
		self.content = json['content']


class Todo:

	def __init__(self, json):
		self.text = json['text']
		self.content = json['content']


class Link:

	def __init__(self, json):
		self.text = json['text']
		self.label = json['label']
		self.type = Type(json['type'])
		try:
			self.signature = json['signature']
		except KeyError:
			self.signature = None
