#!/bin/bash

# User-defined parameters
SAMPLE_NAME=OV606FG
COUNT=10
PROJECT_FILE=/u/laumontc/tsaPaper/${SAMPLE_NAME}/kmerFilt/${SAMPLE_NAME}.${COUNT}.mtechumanAll.0.project
KMER_FILE=/u/laumontc/tsaPaper/${SAMPLE_NAME}/kmerFilt/${SAMPLE_NAME}.${COUNT}.mtechumanAll.0.count
K_LEN=33
ASS_TYPE='linear' # graph --> defines the type of assembly you wanna do
CHECK_RC='FALSE' # TRUE --> check the reverse complement during assembly or not
PATH_OUT=/u/laumontc/tsaPaper/${SAMPLE_NAME}/contigs


# Add required modules
module add boost/1.60.0
module add jellyfish/2.2.3
ulimit -s unlimited


# Create relevant arborescence to store results
ASS_NAME=`echo "${KMER_FILE%.count}" | awk -F [/] '{print $NF}'`

if [ ! -d "$PATH_OUT/" ]; then
    # echo 'Does not exists'
    mkdir $PATH_OUT
    cd $PATH_OUT
    mkdir $ASS_TYPE'_'$ASS_NAME
    cd $ASS_TYPE'_'$ASS_NAME
else
    # echo 'Exists'
    cd $PATH_OUT
    mkdir $ASS_TYPE'_'$ASS_NAME
    cd $ASS_TYPE'_'$ASS_NAME
fi


# only keep kmer and cancer count
# assembly does not like to have too many columns
cut -f2-3 $KMER_FILE > ./$ASS_NAME.assembly


if [ $ASS_TYPE == "linear" ] && [ $CHECK_RC == 'FALSE' ]; then
    echo 'Starting linear assembly - No reverse complement checking';
    /u/laumontc/tsaPaper/scripts/lib/nektar/nektar_cpp/kmer_assembly -p $PROJECT_FILE -k ./$ASS_NAME.assembly \
						-o $PATH_OUT/$ASS_TYPE'_'$ASS_NAME/ -a 1 -l $K_LEN;

elif [ $ASS_TYPE == "linear" ] && [ $CHECK_RC == 'TRUE' ]; then
    echo 'Starting linear assembly - Reverse complement checking';
    /u/laumontc/tsaPaper/scripts/lib/nektar/nektar_cpp/kmer_assembly -r -p $$PROJECT_FILE -k ./$ASS_NAME.assembly \
						-o $PATH_OUT/$ASS_TYPE'_'$ASS_NAME/ -a 1 -l $K_LEN;

elif [ $ASS_TYPE == "shortestPath" ] && [ $CHECK_RC == 'FALSE' ]; then
    echo 'Starting graph assembly - No reverse complement checking';
    /u/laumontc/tsaPaper/scripts/lib/nektar/nektar_cpp/kmer_assembly -p $$PROJECT_FILE -k ./$ASS_NAME.assembly \
						-o $PATH_OUT/$ASS_TYPE'_'$ASS_NAME/ -a 2 -l $K_LEN;

elif [ $ASS_TYPE == "shortestPath" ] && [ $CHECK_RC == 'TRUE' ]; then
    echo 'Starting graph assembly - Reverse complement checking';
    /u/laumontc/tsaPaper/scripts/lib/nektar/nektar_cpp/kmer_assembly -r -p $$PROJECT_FILE -k ./$ASS_NAME.assembly \
						-o $PATH_OUT/$ASS_TYPE'_'$ASS_NAME/ -a 2 -l $K_LEN;

else
    echo 'ASS_TYPE variable must be set to "linear" or "shortestPath"';
    echo 'CHECK_RC variable must be set to "TRUE" or "FALSE"';
fi
