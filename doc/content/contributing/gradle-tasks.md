+++
title = "Gradle Tasks"
weight = 1
+++

## Running Gradle Tasks

Your IDE should give you tools to run Gradle Tasks easily.

For example, in IntelliJ IDEA, there is a Gradle
panel in the upper right of your window. Clicking on that shows all the tasks available. Double-clicking
one of those will run the selected task.

If for some reason, your IDE can't run the Gradle task, you can also use the command line.
In Linux and OSX, open a terminal in the OSPREY project folder and run:
```shell
./gradlew $TASK
```
where `$TASK` is the name of the Gradle task you want to run.

In windows, the command looks slightly different:
```shell
gradlew.bat $TASK
```
