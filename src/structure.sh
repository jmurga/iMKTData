#!/usr/bin/bash

################MANUAL USAGE OPTIONS################
if [ "$1" == "-h" ] || [ "$1" == "-help" ] || [ "$1" == "--h" ] || [ "$1" == "--help" ] || [ -z "$1" ]; then
  printf "\nUsage: `basename $0` execute a simple script to create the main folders structure to start a project.\nPlease to execute it properly just add a name/date to get it as \$1 variable.\n\n\tExample: bash structure.sh 201811\n\n"
  exit 0
fi


mkdir -p $1
mkdir -p $1/scripts $1/rawData $1/results $1/scripts $1/scripts/notebooks
