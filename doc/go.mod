module local.tld/osprey

go 1.17

// NOTE: version numbers have to be in `vx.x.x` format, or go dosen't do the replacement correctly!
require local.tld/hugo-theme-learn v0.0.0
require local.tld/code-python v0.0.0
require local.tld/code-java v0.0.0
require local.tld/code-kotlin v0.0.0

replace local.tld/hugo-theme-learn v0.0.0 => ../build/doc/hugo-theme-learn
replace local.tld/code-python v0.0.0 => ../build/doc/code-python
replace local.tld/code-java v0.0.0 => ../build/doc/code-java
replace local.tld/code-kotlin v0.0.0 => ../build/doc/code-kotlin
