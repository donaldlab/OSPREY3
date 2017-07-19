
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


def init(c):
	wrapStrandFlex(c)
	wrapResidueFlex(c)
	wrapGMECFinder(c)


def wrapStrandFlex(c):
	jtype = jvm.getInnerClass(c.confspace.Strand, 'Flexibility')

	# use array access to call get()
	def getItem(self, key):
		return self.get(key)
	jtype.__getitem__ = getItem


def wrapResidueFlex(c):
	jtype = jvm.getInnerClass(c.confspace.Strand, 'ResidueFlex')

	# wrap setLibraryRotamers() to accept varargs
	def newSetLibraryRotamers(old, self, *mutations):
		return old(self, mutations)
	wrapMethod(jtype, 'setLibraryRotamers', newSetLibraryRotamers)

	# autocast setContinuous() to float
	def newSetContinuous(old, self, angle=None):
		if angle is None:
			return old(self)
		else:
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


def wrapGMECFinder(c):
	jtype = c.gmec.SimpleGMECFinder

	def newFind(old, self, windowSize=None):
		# autocast the find() energy to a float
		# otherwise, jpype will complain if windowSize is really an int
		if windowSize is not None:
			return old(self, float(windowSize))
		else:
			return old(self)
	wrapMethod(jtype, 'find', newFind)

