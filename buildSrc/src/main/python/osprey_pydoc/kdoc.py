
import pathlib
import json


# load the parsed javadoc
with open('build/doc/kdoc/kdoc.json', 'r') as f:
	_kdoc = json.load(f)


default_package = 'edu.duke.cs.osprey'


def _or_empty(str):
	if str is not None:
		return str
	else:
		return ''


class Id:

	def __init__(self, id: str):

		# dokka ID is something like, eg:
		#  package/Class
		#  p1.p2/Class
		#  p1.p2/Class.InnerClass
		#  package/Class/member
		#  package/Class/member/signature
		#  package/Class.Inner1.Inner2/member
		#  .relative.package/Class

		self.raw = id

		# split the path parts
		parts = id.split('/')

		# parse the package, if any, and add the prefix if needed
		if len(parts) >= 1:
			pack = parts[0]
			if pack == '':
				pack = None
			elif pack[0] == '.':
				pack = '%s.%s' % (default_package, pack[1:])
			self.package = pack
		else:
			self.package = None

		# parse the classname
		if len(parts) >= 2:
			cname = parts[1]
			if cname == '':
				cname = None
			self.classname = cname
		else:
			self.classname = None

		# parse the name
		if len(parts) >= 3:
			name = parts[2]
			if name == '':
				name = None
			self.name = name
		else:
			self.name = None

		# parse the signature, if any, and add package prefixes where needed
		if len(parts) >= 4:
			sig_parts = parts[3].split('#')
			for i in range(len(sig_parts)):
				p = sig_parts[i]
				if len(p) > 0 and p[0] == '.':
					p = '%s.%s' % (default_package, p[1:])
				sig_parts[i] = p
			self.signature = '#'.join(sig_parts)
		else:
			self.signature = None


	def __str__(self):
		return 'Id[package=%s, classname=%s, name=%s, signature=%s, raw=%s]' % (self.package, self.classname, self.name, self.signature, self.raw)

	def __repr__(self):
		return self.raw


class Class:

	def __init__(self, id, json):

		self.id = id
		self._json = json

		self.type = Type(json['type'])
		try:
			self.kdoc = Kdoc(json['kdoc'])
		except KeyError:
			self.kdoc = None

		try:
			self.prop_names = json['props']
		except KeyError:
			self.prop_names = []

		try:
			self.func_names = json['funcs']
		except KeyError:
			self.func_names = []

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
		return get_prop('%s/%s' % (self.id.raw, name))


	def prop_or_throw(self, name):
		prop = self.prop(name)
		if prop is None:
			raise Exception('unknown kotlin property: %s, try one of:\n\t%s' % (name, '\n\t'.join(self.prop_names)))
		return prop


	def func(self, name):
		return get_func('%s/%s' % (self.id.raw, name))


	def func_or_throw(self, name):
		return get_func_or_throw('%s/%s' % (self.id.raw, name))


def get_class(id: Id):
	try:
		c = _kdoc['classlikes'][id.package + '/' + id.classname]
	except KeyError:
		return None
	return Class(id, c)


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

	def __init__(self, id, json):

		self.id = id
		self._json = json

		self.name = json['name']
		self.type = Type(json['type'])
		try:
			self.kdoc = Kdoc(json['kdoc'])
		except KeyError:
			self.kdoc = None


def get_prop(id: Id):
	try:
		json = _kdoc['props']['%s/%s/%s/#' % (id.package, _or_empty(id.classname), id.name)]
	except KeyError:
		return None
	return Prop(id, json)


def get_prop_or_throw(id: Id):
	prop = get_prop(id)
	if prop is None:
		raise Exception('unknown kotlin property: %s' % id)
	return prop


class Func:

	def __init__(self, id, json):

		self.id = id
		self._json = json

		self.name = json['name']
		self.args = [FuncArg(a) for a in json['args']]
		self.returns = Type(json['returns'])

		try:
			self.receiver = Type(json['receiver'])
		except KeyError:
			self.receiver = None

		self.url = json['url']

		try:
			self.kdoc = Kdoc(json['kdoc'])
		except KeyError:
			self.kdoc = None


	def find_arg(self, name):
		for arg in self.args:
			if arg.name == name:
				return arg
		return None


	def find_arg_or_throw(self, name):
		arg = self.find_arg(name)
		if arg is None:
			raise Exception('no argument named %s found in function %s\n\ttry one of: %s' % (name, self.name, [a.name for a in self.args]))
		return arg


class FuncArg:

	def __init__(self, json):
		self.name = json['name']
		self.type = Type(json['type'])


def get_func(id: Id):
	try:
		return get_func(id)
	except:
		return None


def get_func_or_throw(id: Id):
	if id.signature is not None:

		# look for overload directly
		try:
			json = _kdoc['funcs']['%s/%s/%s/%s' % (id.package, _or_empty(id.classname), id.name, id.signature)]
		except KeyError:
			raise Exception('unknown kotlin function: %s' % id)

		return Func(id, json)

	else:

		# try to look up the overload
		try:
			overloads = _kdoc['funcs']['%s/%s/%s' % (id.package,  _or_empty(id.classname), id.name)]
		except KeyError:

			if id.classname is None:
				raise Exception('unknown kotlin function: %s' % id)

			# try looking inside the class for suggestions
			c = get_class(id)
			if c is None:
				raise Exception('unknown kotlin class: %s' % id)

			raise Exception('unknown kotlin function `%s` in class `%s`, try one of:\n\t%s' % (id.name, id.classname, '\n\t'.join(c.func_names)))

		if len(overloads) != 1:
			ids = ['%s/%s' % (id.raw, s) for s in overloads]
			raise Exception('kotlin function `%s` has multiple overloads, try one of:\n\t%s' % (id.name, '\n\t'.join(ids)))

		# just one overload, add the signature and try again
		id = Id('%s/%s' % (id.raw, overloads[0]))
		return get_func_or_throw(id)


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
