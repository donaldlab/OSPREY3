
# Osprey

<!--
use raw HTML here, since apparently markdown doesn't have any markup for underlines
-->
<style>
    .u {
        text-decoration: underline;
    }
</style>
<h4>
    <span class="u">O</span>pen-<span class="u">S</span>ource
    <span class="u">P</span>rotein
    <span class="u">Re</span>design
    For
    <span class="u">Y</span>ou
</h4>

OSPREY is developed and maintained by the [Donald Lab][dlab] in the
[Department of Computer Science][dukecs] at [Duke University][duke].

[dlab]: http://www.cs.duke.edu/donaldlab/home.php
[dukecs]: http://www.cs.duke.edu/
[duke]: https://www.duke.edu/

For an introduction to OSPREY 3.0 and its new features please read this paper:

*Journal of Computational Chemistry* 2018; 39(30): 2494-2507.
 * [https://doi.org/10.1002/jcc.25522](https://doi.org/10.1002/jcc.25522) 
 * [Journal Cover Image](http://www.cs.duke.edu/brd/papers/jcc18-osprey3point0/cover-jcc.25043.pdf)
 * [Full-text PDF](http://www.cs.duke.edu/brd/papers/jcc18-osprey3point0/jcc18-osprey-donald.pdf)


## Downloads and Installation

### OSPREY Desktop

OSPREY Desktop is a native desktop application that allows you to pepare designs to run later
on OSPREY server. OSPREY Desktop gives you an interactive graphical interface to:
 * Import molecules from PDB files
 * Edit molecules and add missing information, if any
 * Define conformational flexibility
 * Choose mutations
 * Check for any issues before starting the lengthy design computations

Choose this download to help you prepape your designs for OSPREY Server.

<!-- these <span> HTML tags' contents are auto-generated. if you edit them manually, your changes will be lost -->

#### Linux
* **Download:** <span id="download/desktop/linux/latest">[osprey-desktop-4.0-amd64.deb](https://www.cs.duke.edu/donaldlab/software/osprey/releases/osprey-desktop-4.0-amd64.deb)</span>
* **Installation:** [Installation Instructions for OSPREY Desktop on Linux](install/desktop-linux)

#### OSX
* **Download:** <span id="download/desktop/osx/latest">[osprey-desktop-4.0.dmg](https://www.cs.duke.edu/donaldlab/software/osprey/releases/osprey-desktop-4.0.dmg)</span>
* **Installation:** [Installation Instructions for OSPREY Desktop on OSX](install/desktop-osx)

#### Windows
* **Dowload:** <span id="download/desktop/windows/latest">[osprey-desktop-4.0.msi](https://www.cs.duke.edu/donaldlab/software/osprey/releases/osprey-desktop-4.0.msi)</span>
* **Installation:** [Installation Instructions for OSPREY Desktop on Windows](install/desktop-windows)


[Download older versions](install/versions/#osprey-desktop)


### OSPREY Server

OSPREY Server is a [Python][python] library that provides scripting access to all of OSPREY's features.
OSPREY Server allows you to write Python scripts to:
 * Perform any task available in OSPREY Desktop
 * Run design algorithms on your prepared designs using the best-available computing hardware, including:
   * Multi-core CPUs
   * [CUDA][cuda]-capable GPUs
   * [SLURM][slurm]-capable compute clusters

[python]: https://www.python.org/
[cuda]: http://www.nvidia.com/cuda
[slurm]: https://slurm.schedmd.com/overview.html

Despite the name, OSPREY Server can also be used on desktop and laptop machines.
Although OSPREY is very resource-hungry and works best on powerful machines with lots of compute resources,
particularly CPUs, GPUs, and RAM.

Choose this download if you want to run OSPREY designs and are familiar with Python script programming.

<!-- these <span> HTML tags' contents are auto-generated. if you edit them manually, your changes will be lost -->

#### Linux
 * **Download:** <span id="download/server/linux/latest">[osprey-server-linux-4.0.tbz2](https://www.cs.duke.edu/donaldlab/software/osprey/releases/osprey-server-linux-4.0.tbz2)</span>
 * **Installation:** [Installation Instructions for OSPREY Server on Linux](install/server-linux)

#### OSX
 * **Download:** <span id="download/server/osx/latest">[osprey-server-osx-4.0.tbz2](https://www.cs.duke.edu/donaldlab/software/osprey/releases/osprey-server-osx-4.0.tbz2)</span>
 * **Installation:** [Installation Instructions for OSPREY Server on OSX](install/server-osx)

#### Windows
 * **Dowload:** <span id="download/server/windows/latest">[osprey-server-win-4.0.tbz2](https://www.cs.duke.edu/donaldlab/software/osprey/releases/osprey-server-win-4.0.tbz2)</span>
 * **Installation:** [Installation Instructions for OSPREY Server on Windows](install/server-windows)


[Download older versions](install/versions/#osprey-server)


### OSPREY Service

Although OSPREY is a cross-platform application, not all of the design preparation features work on all platforms.
OSPREY Service is an internet service that provides these features for platforms where they are not supported.
The Donald Lab provides a default service that OSPREY will use automatically when preparing designs, so
most people will not need to download OSPREY service. Hovewer, if you don't want to rely on the service
provided by the Donald Lab, then you can self-host your own service.

The OSPREY Service is only supported on Linux.

Choose this download to self-host your own service for OSPREY design preparation.

<!-- these <span> HTML tags' contents are auto-generated. if you edit them manually, your changes will be lost -->

#### Linux
* **Download:** <span id="download/service-docker/linux/latest">[osprey-service-docker-0.3.tar](https://www.cs.duke.edu/donaldlab/software/osprey/releases/osprey-service-docker-0.3.tar)</span>
* **Installation:** [Installation Instructions for OSPREY Service on Linux](install/service-docker-linux)


[Download older versions](install/versions/#osprey-service)


## Citing OSPREY

If you publish or present results from OSPREY, please mention the name "OSPREY"
and cite our papers as described in [CITING_OSPREY.txt][citing].

[citing]: https://github.com/donaldlab/OSPREY3/blob/main/CITING_OSPREY.txt


## License

GPLv2

Copyright (C) 2017 Duke University

This program is free software; you can redistribute it and/or modify it under the terms of the
GNU General Public License version 2 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

The full text of the GPLv2 is included in [LICENSE.txt][license].

[license]: https://github.com/donaldlab/OSPREY3/blob/main/LICENSE.txt
