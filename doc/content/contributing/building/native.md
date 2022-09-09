
+++
menuTitle = "Native Code"
title = "Building OSPREY's Native Code"
weight = 4
+++


This section explains how to build some OSPREY modules that are written in
"native" languages like C++ and its variants like CUDA.
In this context, "native" means code that runs directly on the operating system
and not the Java Virtual Machine (JVM).


## Native libaries

To make life easier for Java developers, the build artifacts for OSPREY's native libraries
are saved in the git repository. This means developers who want to build OSPREY usually only
need to use the usual Gradle tasks to make builds, and they don't have to worry about build
toolchains for other languages or ecosystems.

You would only need to explicitly build a native module if you were making changes to the
native code itself. However, if you do make changes to a native module, make sure the final
build artifact is built in release mode and committed to the git repository when you're done.


## Using CMake

Most of OSPREY's native code is built using [CMake][cmake].

[cmake]: https://cmake.org/

To build CMake projects, follow the usual incantation for CMake:

* Make a build folder to do an out-of-source build. eg, `cmake-build-debug` or `cmake-build-release`.
  {{% notice note %}}
  The `.gitignore` files for OSPREY native modules are set to ignore folders matching the glob `cmake-*`.
  {{% /notice %}}
* `cd` to the build folder in your favorite terminal.
* run `cmake ..` to do some first-time setup.
* run `make $TARGET` to do a build, where `$TARGET` is the make target.

CMake builds are configured in the `CMakeLists.txt` file at the module root folder.


## Native Modules in OSPREY

### Conformation Energy Calculator

**Location:** `/src/main/cc/ConfEcalc`

Useful CMake targets:

 * **`ConfEcalc`:**
   This task builds the native library.

 * **`ConfEcalc_CopyLibs`:**
   This task calls `ConfEcalc` and then copies the built library to `/src/main/resources/$ARCH`
   where the [JNA][jna] loader in the Java code expects to find it.
   `$ARCH` is a unique string for your local platform, eg `linux-x86_64`, `darwin--aarch64`.
   Use this task if you want to run the Java code using your newly-built library.

[jna]: https://github.com/java-native-access/jna


### CUDA Conformation Energy Calculator

**Location:** `/src/main/cu/CudaConfEcalc`

Useful CMake targets:

* **`CudaConfEcalc`:**
  This task builds the native library and the CUDA kernels.

* **`CudaConfEcalc_CopyLibs`:**
  This task calls `CudaConfEcalc` and then copies the built library to `/src/main/resources/$ARCH`
  where the [JNA][jna] loader in the Java code expects to find it.
  `$ARCH` is a unique string for your local platform, eg `linux-x86_64`, `darwin--aarch64`.
  Use this task if you want to run the Java code using your newly-built library.
