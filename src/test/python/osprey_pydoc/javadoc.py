
import pathlib
import json


# load the parsed javadoc
with open('build/javadoc.json', 'r') as f:
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

	def __init__(self, path, json):
		self.type = Type(json['type'])
		try:
			self.javadoc = json['javadoc']
		except KeyError:
			self.javadoc = None


def get_class(path: Path):
	try:
		c = _javadoc[path.classname]
	except KeyError:
		return None
	return Class(path, c)


class Field:

	def __init__(self, path, json):
		self.name = path.member
		self.type = Type(json['type'])
		try:
			self.javadoc = json['javadoc']
		except KeyError:
			self.javadoc = ''
		try:
			self.initializer = json['initializer']
		except KeyError:
			self.initializer = None


def get_field(path: Path):
	try:
		field = _javadoc[path.classname]['fields'][path.member]
	except KeyError:
		return None
	return Field(path, field)


class Method:

	def __init__(self, path, json):
		self.name = path.member
		self.signature = json['signature']
		try:
			self.returns = Type(json['returns'])
		except KeyError:
			self.returns = None
		try:
			self.javadoc = json['javadoc']
		except KeyError:
			self.javadoc = None


def get_method(path: Path):
	try:
		method = _javadoc[path.classname]['methods'][path.member]
	except KeyError:
		return None
	return Method(path, method)


class Type:

	def __init__(self, json):
		try:
			self.name = json['name']
		except (TypeError, KeyError):
			# no object here, the json bit must be the string name
			self.name = json
			return
		try:
			self.url = json['url']
		except KeyError:
			self.url = None
		try:
			self.params = [Type(param) for param in json['params']]
		except KeyError:
			self.params = None
