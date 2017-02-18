
import os
import docutils, sphinx
import javalang


ast_cache = {}


# Sphinx application API
# www.sphinx-doc.org/en/1.5.1/extdev/appapi.html

def setup(app):

	# configs
	app.add_config_value('javadoc_sources_dir', '.', 'env')
	app.add_config_value('javadoc_package_prefix', '', 'env')

	app.add_domain(JavadocDomain)

	app.add_role('java:fielddoc', role_fielddoc)

	return { 'version': '0.1' }


def expand_classname(classname, config):
	
	# expand the classname if we need
	if classname[0] == '.':
		classname = config.javadoc_package_prefix + classname

	return classname


def read_file(path):
	with open(path, 'r') as file:
		return file.read()


def get_class_ast(classname, config):

	# get the outer class
	outer_classname = classname.split('$')[0]

	# convert to a filepath
	path = outer_classname.replace('.', '/')
	path = '%s/%s.java' % (config.javadoc_sources_dir, path)

	# read from the cache, or read from the file
	try:
		return ast_cache[path]
	except KeyError:

		ast = javalang.parse.parse(read_file(path))
		ast_cache[path] = ast

		return ast
	

def strip_javadoc(javadoc):

	if javadoc[0:3] == '/**':
		javadoc = javadoc[3:]
	if javadoc[-2:] == '*/':
		javadoc = javadoc[:-2]
	return javadoc.strip()


def is_public(thing):
	return 'public' in thing.modifiers


def make_method_signature(method):
	return method.name


class ParsingDirective(docutils.parsers.rst.Directive):

	@property
	def env(self):
		return self.state.document.settings.env


	@property
	def config(self):
		return self.env.config


	def parse(self, rst):

		# if rst is a list of strings, join it
		if not isinstance(rst, str):
			rst = '\n'.join(rst)

		# parse the rst and return the nodes
		docSource = 'I come from the water'
		doc = docutils.utils.new_document(docSource, self.state.document.settings)
		docutils.parsers.rst.Parser().parse(rst, doc)
		return doc.children


class JavaClassDirective(ParsingDirective):
	
	has_content = True

	def get_class_ast(self, classname):
		classname = expand_classname(classname, self.config)
		return get_class_ast(classname, self.config)


	def run(self):

		ast = self.get_class_ast(self.content[0])

		# use the top class
		classtype = ast.types[0]

		rst = []

		# class name header
		rst.append(classtype.name)
		rst.append('=' * len(classtype.name))

		rst.append('')

		# show fields
		rst.append('Attributes')
		rst.append('----------')
		rst.append('')
		for field in classtype.fields:
			
			# skip fields without javadoc or aren't public
			if field.documentation is None or not is_public(field):
				continue

			# show the field name and javadoc
			for decl in field.declarators:
				rst.append('.. py:attribute:: %s' % decl.name)
				rst.append('')
				rst.append('\t' + strip_javadoc(field.documentation))
				rst.append('')

		# show methods
		rst.append('Methods')
		rst.append('-------')
		rst.append('')
		for method in classtype.methods:
			
			# skip methods without javadoc or aren't public
			if method.documentation is None or not is_public(method):
				continue

			rst.append('.. py:method:: %s' % make_method_signature(method))
			rst.append('')
			rst.append('\t' + strip_javadoc(method.documentation))
			rst.append('')

		return self.parse(rst)


class JavaRole(sphinx.roles.XRefRole):

	def process_link(self, env, refnode, has_explicit_title, title, target):

		# nothing to do here, all the work can be done in Domain.resolve_xref()
		return title, target


# how to write function roles:
# http://docutils.sourceforge.net/docs/howto/rst-roles.html

def role_fielddoc(name, rawtext, text, lineno, inliner, options={}, content=[]):
	
	env = inliner.document.settings.env

	# split the classname and the field name
	text_parts = text.split('.')
	classname = '.'.join(text_parts[:-1])
	full_classname = expand_classname(classname, env.config)
	fieldname = text_parts[-1]

	# find the field javadoc
	ast = get_class_ast(full_classname, env.config)
	def find_javadoc():
		for field in ast.types[0].fields:
			for declarator in field.declarators:
				if declarator.name == fieldname:
					return field.documentation
		return '(no found javadoc for field: %s)' % text
	javadoc = find_javadoc()

	# make the text node
	javadoc = strip_javadoc(javadoc)
	return [docutils.nodes.Text(javadoc)], []
	

# Sphinx domain API:
# http://www.sphinx-doc.org/en/1.5.1/extdev/domainapi.html

class JavadocDomain(sphinx.domains.Domain):

	name = 'java'
	label = 'Osprey javadoc processor'

	directives = {
		'class': JavaClassDirective
	}

	roles = {
		'class': JavaRole()
	}

	def resolve_xref(self, env, fromdocname, builder, typ, target, node, contnode):
		
		# get the fully-qualified java class name and the short class name
		full_classname = expand_classname(target, env.config)
		short_classname = full_classname.split('.')[-1]

		# reset the next of the link
		contnode.replace(
			contnode.children[0],
			docutils.nodes.Text(short_classname)
		)

		# remove any leading . from the target
		if target[0] == '.':
			target = target[1:]

		# build the ref node
		path = 'api.%s' % target
		anchor = ''
		title = None

		return sphinx.util.nodes.make_refnode(
			builder, fromdocname, path, anchor, contnode, title
		)
 
