# this is the scripts used to do RAA-seq rutine analysis from raw sequencing data


#!/usr/bin/bash
################################################################################
# Help                                                                         #
################################################################################
if [ "$1" == "-h" | "$1" == "--help"] 
then
    echo "Usage 1: ./ `basename $0`  <anyname>"
    echo "Usage 2: sh `basename $0` -a <anyname>"
    echo "This Help File: sh `basename $0` -h"
    exit 0
fi



