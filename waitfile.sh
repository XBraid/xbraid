#!/bin/bash
file=$1
sz=$(stat --printf="%s" ${file})
nsz=${sz}
while [[ "${sz}" == "${nsz}" ]]
do
    sleep 0.1
    nsz=$(stat --printf="%s" ${file})
done
[[ "${sz}" != "${nsz}" ]] && echo ${file}
