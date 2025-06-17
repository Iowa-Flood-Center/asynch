#!/bin/sh

# This script should not be needed, but apparently docker needs it to run quotes inside commands, so well, it is only one instruction.

cd ../build
../configure CFLAGS="-O3 -DNDEBUG -Wno-format-security" --disable-dependency-tracking


