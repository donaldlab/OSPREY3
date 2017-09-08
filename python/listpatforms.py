
import pip.pep425tags

for (pyver, abi, platform) in pip.pep425tags.supported_tags:
	print('%s.%s.%s' % (pyver, abi, platform))

