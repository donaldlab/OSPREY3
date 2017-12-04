
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
	wrapStrandBuilder(c)


def checkDeprecatedResNumbers(resNums, fnnames):

	import traceback

	hasIntResNum = False
	for key in resNums:
		if isinstance(key, (int, long, float)):
			hasIntResNum = True
			break

	if hasIntResNum:

		print('WARNING: Number-valued residue numbers (e.g., 4, 123, 42) are deprecated and will be not be supported in a future version.'
			+ ' Instead, use String-valued residue numbers prefixed with the chain ID if available (e.g., \'A4\', \'G123\', \'42\').')

		# get the stack trace, and remove the specified osprey frames
		stacktrace = traceback.extract_stack()
		for i in range(len(fnnames)):
			fnname = stacktrace[-1][2]
			if fnname == fnnames[i]:
				stacktrace = stacktrace[:-1]

		print('\tresidue numbers: %s' % resNums)
		for (filename, lineno, fnname, line) in stacktrace:
			print('\tat %s:%d' % (filename, lineno))
			#print((filename, lineno, fnname))


def wrapStrandFlex(c):
	jtype = jvm.getInnerClass(c.confspace.Strand, 'Flexibility')

	# use array access to call get()
	def getItem(self, key):

		# TEMP: check for integer-valued residue numbers, which are now deprecated
		checkDeprecatedResNumbers([key], ['checkDeprecatedResNumbers', 'getItem'])

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


def wrapStrandBuilder(c):
	jtype = jvm.getInnerClass(c.confspace.Strand, 'Builder')

	# public Builder setResidues(int firstResNum, int lastResNum) {
	def newSetResidues(old, self, start, stop):

		# TEMP: check for integer-valued residue numbers, which are now deprecated
		checkDeprecatedResNumbers([start, stop], ['checkDeprecatedResNumbers', 'newSetResidues', 'curried', 'Strand'])

		old(self, start, stop)
	wrapMethod(jtype, 'setResidues', newSetResidues)
