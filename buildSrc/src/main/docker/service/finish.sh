#!/bin/sh

# this script is called when an s6 service exist

# when the service dies, tell s6 to exit and not to restart it
# ideally, the whole container should die if any one service dies
kill -INT 1

# tell s6 not to restart this service
exit 125
