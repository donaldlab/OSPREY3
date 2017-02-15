
import jvm


def wrapMethod(obj, name, newMethod):
	oldMethod = getattr(obj, name)
	def curried(self, *args):
		return newMethod(oldMethod, self, *args)
	setattr(obj, name, curried)


def wrapProperty(obj, name, getter, setter=None):

	oldProp = getattr(obj, name)

	def tempGetter(self):
		return getter(oldProp.__get__, self)
	def tempSetter(self, oldProp, val):
		setter(oldProp.__set__, self, val)

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
	wrapStrandFlex()
	wrapResidueFlex()
	wrapGMECFinder()


def wrapStrandFlex():
	jtype = getJavaClass('confspace.Strand$Flexibility')

	# use array access to call get()
	def getItem(self, key):
		return self.get(key)
	jtype.__getitem__ = getItem


def wrapResidueFlex():
	jtype = getJavaClass('confspace.Strand$ResidueFlex')

	# wrap setLibraryRotamers() to accept varargs
	def newSetLibraryRotamers(old, self, *mutations):
		return old(self, jvm.toArrayList(mutations))
	wrapMethod(jtype, 'setLibraryRotamers', newSetLibraryRotamers)

	# autocast setContinuous() to float
	def newSetContinuous(old, self, angle):
		return old(self, float(angle))
	wrapMethod(jtype, 'setContinuous', newSetContinuous)

	# autocast setContinuousEllipses() to float
	def newSetContinuousEllipses(old, self, angle):
		return old(self, float(angle))
	wrapMethod(jtype, 'setContinuousEllipses', newSetContinuousEllipses)
	
	# NOTE: property wrapping example
	#def newGetter(oldGetter, self):
	#	return oldGetter(self)
	#wrapProperty(strand, 'flexibility', newGetter)


def wrapGMECFinder():
	jtype = getJavaClass('gmec.SimpleGMECFinder')

	def newFind(old, self, windowSize=None):
		# autocast the find() energy to a float
		# otherwise, jpype will complain if windowSize is really an int
		if windowSize is not None:
			return old(self, float(windowSize))
		else:
			return old(self)
	wrapMethod(jtype, 'find', newFind)

