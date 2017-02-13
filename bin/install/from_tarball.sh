#!/usr/bin/env bash

mkdir -p src/

pushd src/

wget \
    --mirror \
    --no-host-directories \
    --cut-dirs=2 \
    ftp://ftp.pantherdb.org/hmm_scoring/current_release/pantherScore2.0/

popd