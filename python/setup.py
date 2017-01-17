
import setuptools

setuptools.setup(
	name='osprey',
	version='3.0-alpha1',
	description='Open-Source Protein Redesign for You',
	url='https://github.com/donaldlab/OSPREY_refactor',
	packages=setuptools.find_packages(exclude=['contrib', 'docs', 'tests']),
	install_requires=['JPype1'],
	# TODO: how to package the jar file?
	#package_data={
	#	'sample': 'osprey.jar'
	#},
	classifiers=[
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3',
		'Operating System :: OS Independent'
	]
)

# TODO: get a pypi account and upload our python package when it's ready for release
# then installs can be as easy as `pip install osprey`

