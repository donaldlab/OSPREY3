#!/bin/bash

for master in PM-*/
do

pmk="STILL_RUNNING"
mutk="STILL_RUNNING"

for dir in ${master}*
do

lname="${dir#*\/}"

if [ ! -f "$dir/submit.out" ]; then
    continue
fi

while IFS=',' read -r f1 f2 f3 f4
do

# get PM
if [ "$lname" == "PM" ] && [ "$f1" == " 1" ]
then
pmk="${f3%% *}"

# get mutuant/ortholog
elif [ "$lname" != "PM" ] && [ "$f1" == " 1" ]
then
mutk="${f3%% *}"
mname="$lname"

fi

done < "$dir/submit.out"
done

# still get mut name if we don't have the score
if [ "$mutk" == "STILL_RUNNING" ]
then
holder="${dir%\/*}"
mname="${holder#*PM-}"
fi

echo "$mname,$pmk,$mutk"

done


