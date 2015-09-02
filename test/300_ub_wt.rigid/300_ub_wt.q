#!/bin/csh -f
# ---------------------------
# our name 
#$ -N test300_ub_wt
#
#
# Operate in current working directory
#$ -cwd
#$ -S /bin/csh
#
#$ -l mem_free=16G
#$ -o 300_ub_wt.out -j y
#
# ---------------------------
set jdep=../lib
set jars=$jdep/joptimizer.jar:$jdep/architecture-rules-3.0.0-M1.jar:$jdep/colt-1.2.0.jar:$jdep/commons-beanutils-1.6.jar:$jdep/commons-collections-2.1.jar:$jdep/commons-digester-1.6.jar:$jdep/commons-io-1.4.jar:$jdep/commons-lang-2.5.jar:$jdep/commons-logging-1.1.1.jar:$jdep/commons-math3-3.0.jar:$jdep/jdepend-2.9.1.jar:$jdep/junit-3.8.1.jar:$jdep/log4j-1.2.14.jar:$jdep/xml-apis-1.0.b2.jar

java -Xmx16g -cp build/classes:$jars edu.duke.cs.osprey.control.Main -c KStar.cfg findGMEC System_300_ub.cfg DEE_300_ub_wt.cfg
