
import pathlib
import json


# load the parsed javadoc
with open('build/doc/kdoc.json', 'r') as f:
	_kdoc = json.load(f)


default_package = 'edu.duke.cs.osprey'


class Id:

	def __init__(self, id: str):

		# dokka ID is something like, eg:
		#  package/Class
		#  p1.p2/Class
		#  p1.p2/Class.InnerClass
		#  package/Class/member
		#  package/Class.Inner1.Inner2/member
		#  .relative.package/Class

		# apply the default package root if needed
		if id[0] == '.':
			id = '%s.%s' % (default_package, id[1:])

		self.raw = id

		# split the path parts
		parts = id.split('/')
		if len(parts) >= 2:
			self.package = parts[0]
			self.classname = parts[1]
		else:
			raise Exception('malformed kotlin ID: %s' % id)
		if len(parts) > 2:
			self.name = parts[2]
		else:
			self.name = None


	def __str__(self):
		return 'Id[package=%s, classname=%s, name=%s]' % (self.package, self.classname, self.name)

	def __repr__(self):
		return self.raw


class Class:

	def __init__(self, json):

		self._json = json

		self.type = Type(json['type'])
		try:
			self.kdoc = Kdoc(json['kdoc'])
		except KeyError:
			self.kdoc = None

		try:
			self.prop_names = [p['name'] for p in json['props']]
		except KeyError:
			self.prop_names = []

		# TODO: functions

		try:
			self.entry_names = [f['name'] for f in json['entries']]
		except KeyError:
			self.entry_names = []


	def entry(self, name):
		try:
			i = self._json['entriesLut'][name]
		except KeyError:
			return None
		return Entry(self._json['entries'][i])


	def entry_or_throw(self, name):
		entry = self.entry(name)
		if entry is None:
			raise Exception('unknown kotlin enum entry: %s, try one of:\n\t%s' % (name, '\n\t'.join(self.entry_names)))
		return entry


	def prop(self, name):
		try:
			i = self._json['propsLut'][name]
		except KeyError:
			return None
		return Prop(self._json['props'][i])


	def prop_or_throw(self, name):
		prop = self.prop(name)
		if prop is None:
			raise Exception('unknown kotlin property: %s, try one of:\n\t%s' % (name, '\n\t'.join(self.prop_names)))
		return prop


def get_class(id: Id):
	try:
		c = _kdoc[id.package + '/' + id.classname]
	except KeyError:
		return None
	return Class(c)


def get_class_or_throw(id: Id):
	c = get_class(id)
	if c is None:
		raise Exception('unknown kotlin class: %s' % id)
	return c


class Entry:

	def __init__(self, json):

		self._json = json

		self.name = json['name']
		try:
			self.kdoc = Kdoc(json['kdoc'])
		except KeyError:
			self.kdoc = None


class Prop:

	def __init__(self, json):

		self._json = json

		self.name = json['name']
		self.type = Type(json['type'])
		try:
			self.kdoc = Kdoc(json['kdoc'])
		except KeyError:
			self.kdoc = None


def get_prop(id: Id):
	if id.classname is None:

		# look in the top-level
		try:
			json = _kdoc[id.package + '//' + id.name]
		except KeyError:
			return None
		return Prop(json)

	else:

		# look in the enclosing class
		c = get_class(id)
		if c is None:
			return None
		return c.prop(id.name)


def get_prop_or_throw(id: Id):
	if id.classname is None:

		# look in the top-level
		try:
			json = _kdoc[id.package + '//' + id.name]
		except KeyError:
			raise Exception('unknown kotlin property: %s' % id)
		return Prop(json)

	else:

		# look in the enclosing class
		c = get_class_or_throw(id)
		return c.prop_or_throw(id.name)


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
			self.nullable = json['nullable']
		except (KeyError, TypeError):
			self.nullable = False

		# TODO: do we need to bother with rendering functional types correctly?
		#  do they even appear in the Python API at all?
		try:
			self.functional = json['functional']
		except (KeyError, TypeError):
			self.functional = False

		try:
			self.variance = json['variance']
		except (KeyError, TypeError):
			self.variance = None

		try:
			self.params = [Type(param) for param in json['params']]
		except (KeyError, TypeError):
			self.params = None


class Kdoc:

	def __init__(self, text):
		self.text = text
