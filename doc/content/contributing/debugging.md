+++
title = "Debugging"
weight = 2
+++

## Debugging OSPREY JVM code from Python scripts

While JPype obfuscates Java from the user, Java is front-and-center for the developer.  When bugs or
implementation issues arise, The program will probably need to be debugged in Java. Yet, since
typical input is a Python script, the process of debugging is not as straightforward as it is for a
normal java application. That said, it's also not too complex. The way to do it is to pass an
addition argument to the `osprey.start()` function:

```python
osprey.start(attachJvmDebugger=True)
```

This will enable the Java process to be attached to by a java debugger. [Intellij IDEA][idea] has a nice one you
can use. When you pass this optional parameter to the `start()` function, the script will pause
after the JVM starts and print a message telling you to attach a debugger. Once you hit enter, the
script will resume. If you have not attached a debugger or have attached a debugger but have not
set any breakpoints, the script runs without pause, as before. The process for attaching a debugger
in Intellij is as follows:

[idea]: https://www.jetbrains.com/idea/

1. Open the source code for Osprey as a java project in Intellij. The java and python source code
should be the same version, i.e. it's best to uninstall any existing Osprey installations,
download the Osprey repository's source code, and execute ``./gradlew pythonDevelop`` so that the
python version of Osprey is the same as the java code you're about to debug.
2. Run the Osprey python script with the ``attachJvmDebugger`` argument set to True.
3. The python interpreter will pause, giving you an opportunity to attach the java debugger. In
Intellij, this can be done with the *Attach to Process* option in the *Run* menu.
4. Intellij intelligently only shows java processes. On my computer, the correct one to select is
the one entitled, counterintuitively, ``python``.
5. Set a breakpoint somewhere in code you know will be hit.
6. Return the to the python script and click enter to continue.

Your breakpoint should be hit and you should be on your way to debugging as normal in java land.
