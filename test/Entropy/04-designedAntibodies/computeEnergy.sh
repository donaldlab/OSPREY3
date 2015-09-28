#!/bin/bash
export CLASSPATH=/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/bin
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/colt-1.2.0.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/architecture-rules-3.0.0-M1.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/commons-beanutils-1.6.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/commons-collections-2.1.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/commons-digester-1.6.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/commons-io-1.4.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/commons-lang-2.5.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/commons-logging-1.1.1.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/commons-math3-3.0.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/jdepend-2.9.1.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/joptimizer.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/junit-3.8.1.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/log4j-1.2.14.jar
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/pablo/devel/osprey-GIT/OSPREY/lib/xml-apis-1.0.b2.jar
export CLASSPATH=$CLASSPATH:/home/home1/pablo/dlab/soft/gurobi563/linux64/lib/gurobi.jar
java  -Xmx32000M KStar -c ../KStar.cfg -t 44 computeEMolHIV System.cfg DEE.cfg $1 | tee interfaceEnergy
#java  -Xmx2000M KStar -t 15 genStructDEE System.cfg GenStructDEE.cfg
