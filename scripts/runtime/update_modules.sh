#!/usr/bin/env bash

cwd=${pwd}
cd innovation-lab
git pull --recurse-submodules && cd dodo-cloning-kit && git pull https://github.com/ndbrown6/dodo-cloning-kit.git master
cd ${cwd}
