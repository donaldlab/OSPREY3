
import jvm


def wrapMethod(ctype, name, newMethod):
	oldMethod = getattr(ctype, name)
	def curried(self, *args):
		return newMethod(self, oldMethod, *args)
	setattr(ctype, name, curried)


def wrapProperty(ctype, name, getter, setter=None):

	oldProp = getattr(ctype, name)

	def tempGetter(self):
		return getter(self, oldProp.__get__)
	def tempSetter(self, oldProp, val):
		setter(self, oldProp.__set__, val)

	newProp = None
	if setter is None:
		newProp = property(tempGetter)
	else:
		newProp = property(tempGetter, tempSetter)

	setattr(ctype, name, newProp)


def getClass(classname):
	classname = 'edu.duke.cs.osprey.' + classname
	jclass = jvm.c.java.lang.Class
	classloader = jvm.c.java.lang.ClassLoader.getSystemClassLoader()
	return jclass.forName(classname, True, classloader)
	

def init():
	augmentStrandFlex(getClass('confspace.Strand$Flexibility'))
	augmentResidueFlex(getClass('confspace.Strand$ResidueFlex'))


def augmentStrandFlex(ctype):

	# use array access to call get()
	def getItem(self, key):
		return self.get(key)
	ctype.__getitem__ = getItem


def augmentResidueFlex(ctype):

	# augment setLibraryRotamers() to accept varargs
	def newSetLibraryRotamers(self, old, *mutations):
		old(self, jvm.toArrayList(mutations))
	wrapMethod(ctype, 'setLibraryRotamers', newSetLibraryRotamers)
	
	# NOTE: property wrapping example
	#def newGetter(self, oldGetter):
	#	return oldGetter(self)
	#wrapProperty(strand, 'flexibility', newGetter)

