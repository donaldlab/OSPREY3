
import jvm


def wrapMethod(obj, name, newMethod):
	oldMethod = getattr(obj, name)
	def curried(self, *args):
		return newMethod(self, oldMethod, *args)
	setattr(obj, name, curried)


def wrapProperty(obj, name, getter, setter=None):

	oldProp = getattr(obj, name)

	def tempGetter(self):
		return getter(self, oldProp.__get__)
	def tempSetter(self, oldProp, val):
		setter(self, oldProp.__set__, val)

	newProp = None
	if setter is None:
		newProp = property(tempGetter)
	else:
		newProp = property(tempGetter, tempSetter)

	setattr(obj, name, newProp)


def getJavaClass(classname):
	classname = 'edu.duke.cs.osprey.' + classname
	jclass = jvm.c.java.lang.Class
	classloader = jvm.c.java.lang.ClassLoader.getSystemClassLoader()
	return jclass.forName(classname, True, classloader)
	

def init():
	augmentStrandFlex()
	augmentResidueFlex()
	augmentGMECFinder()


def augmentStrandFlex():
	jtype = getJavaClass('confspace.Strand$Flexibility')

	# use array access to call get()
	def getItem(self, key):
		return self.get(key)
	jtype.__getitem__ = getItem


def augmentResidueFlex():
	jtype = getJavaClass('confspace.Strand$ResidueFlex')

	# augment setLibraryRotamers() to accept varargs
	def newSetLibraryRotamers(self, old, *mutations):
		return old(self, jvm.toArrayList(mutations))
	wrapMethod(jtype, 'setLibraryRotamers', newSetLibraryRotamers)

	# autocast setContinuous() to float
	def newSetContinuous(self, old, angle):
		return old(self, float(angle))
	wrapMethod(jtype, 'setContinuous', newSetContinuous)

	# autocast setContinuousEllipses() to float
	def newSetContinuousEllipses(self, old, angle):
		return old(self, float(angle))
	wrapMethod(jtype, 'setContinuousEllipses', newSetContinuousEllipses)
	
	# NOTE: property wrapping example
	#def newGetter(self, oldGetter):
	#	return oldGetter(self)
	#wrapProperty(strand, 'flexibility', newGetter)


def augmentGMECFinder():
	jtype = getJavaClass('gmec.SimpleGMECFinder')

	def newFind(self, old, windowSize=None):
		# autocast the find() energy to a float
		# otherwise, jpype will complain if windowSize is really an int
		if windowSize is not None:
			return old(self, float(windowSize))
		else:
			return old(self)
	wrapMethod(jtype, 'find', newFind)

