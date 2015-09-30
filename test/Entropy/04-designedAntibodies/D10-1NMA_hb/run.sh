#!/bin/bash


CLASSPATH=/usr/project/dlab/Users/hunter/OSPREY_CODE/OSPREY_refactor_v2/lib/
CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/hunter/OSPREY_CODE/OSPREY_refactor_v2/build/classes/edu/duke/cs/osprey/control/
export CLASSPATH=$CLASSPATH:/usr/project/dlab/Users/hunter/OSPREY_CODE/OSPREY_refactor_v2/dist/OSPREY_refactor_v2.jar
echo $CLASSPATH

java -cp /usr/project/dlab/Users/hunter/OSPREY_CODE/OSPREY_refactor_v2/dist/OSPREY_refactor_v2.jar edu.duke.cs.osprey.control.Main  -c KStar.cfg findGMEC DEE.cfg System.cfg
