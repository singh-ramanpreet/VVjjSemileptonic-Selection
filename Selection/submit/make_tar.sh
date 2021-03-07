#!/bin/bash

cmssw=$(basename $CMSSW_BASE)

rm -f ${cmssw}.tgz

tar --exclude-vcs -zcf ${cmssw}.tgz -C $CMSSW_BASE/.. \
    ${cmssw}/lib \
    ${cmssw}/bin \
    ${cmssw}/python \
    ${cmssw}/config \
    ${cmssw}/.SCRAM \
    ${cmssw}/src/VVjjSemileptonic/Selection/data \
    $(find ${CMSSW_BASE}/src/VVjjSemileptonic/Selection -type f -name *.txt | sed -E 's|.*(CMSSW_)|\1|p')
