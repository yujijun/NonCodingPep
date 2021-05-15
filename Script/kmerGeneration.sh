SAMPLE_NAME=OV642_CD45-Ascites
PROJECT_HASH=d9b605a68f5437045dfc061fa3d9e579
PATH_TO_FASTQ=https://bioinfo.iric.ca/seq/${PROJECT_HASH}
REP=(HL2FJBGX5/fastq/Sample_OV642_CD45-Ascites) #('_102015/' '_wtf')
NB_THREADS=8
# trimmomatic parameters
qLead=20
qTrail=20
minLen=33
# decided to use the longest k-mer size to be used
# avoid including more information in databases generated with shorter ksize

# jellyfish parameters
STRANDED='TRUE' # or FALSE
LEN=(24 33)
HASH_SIZE='1G'
OUTPATH=/u/laumontc/tsaPaper/${SAMPLE_NAME}/genomicData

# # Load required modules
module add jellyfish/2.2.3
module add fastx-toolkit/0.0.14

# # Find max value in length
# max=0
# for i in ${LEN[@]}
# do
# if [[ $i -gt $max ]] ; then
# max=$i
# fi
# done


# Create appropriate file hierarchy
cd $TMPDIR

mkdir $SAMPLE_NAME
cd $SAMPLE_NAME


# # Download fastq files files
echo `date "+%Y-%m-%d %H:%M:%S":` 'Downloading fastq.gz files...'
for r in "${REP[@]}"
do
    echo '    Processing replicate' ${PATH_TO_FASTQ}/${r}'...' `date "+%Y-%m-%d %H:%M:%S"`
    wget -r -l 1 -A 'fastq.gz' -nH --cut-dirs 2 --no-parent -e robots=off ${PATH_TO_FASTQ}/${r}
    # -r --no-parent -nd *.fastq.gz ${PATH_TO_FASTQ}${r} ## Old params
    # -P <path> to specify output directory, by default .
done
echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'
echo ''


# Concatenate files
# R1
echo `date "+%Y-%m-%d %H:%M:%S":` 'Concatenating R1 fastqs...'
\ls -1 ./*/fastq/Sample_${SAMPLE_NAME}/*_R1_*.fastq.gz
cat ./*/fastq/Sample_${SAMPLE_NAME}/*_R1_*.fastq.gz > ./${SAMPLE_NAME}_R1.fastq.gz
rm -v ./*/fastq/Sample_${SAMPLE_NAME}/*_R1_*.fastq.gz
cp ./${SAMPLE_NAME}_R1.fastq.gz ${OUTPATH}
echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'
echo ''
# R2
echo `date "+%Y-%m-%d %H:%M:%S":` 'Concatenating R2 fastqs...'
\ls -1 ./*/fastq/Sample_${SAMPLE_NAME}/*_R2_*.fastq.gz
cat ./*/fastq/Sample_${SAMPLE_NAME}/*_R2_*.fastq.gz > ./${SAMPLE_NAME}_R2.fastq.gz
rm -v ./*/fastq/Sample_${SAMPLE_NAME}/*_R2_*.fastq.gz
cp ./${SAMPLE_NAME}_R2.fastq.gz ${OUTPATH}
echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'
echo ''


# Remove adapters from fastqs
# workinf with TrueSeq adapters
echo `date "+%Y-%m-%d %H:%M:%S":` 'Removing adapters...'
echo '    Parameters used for Trimmomatic:'
echo '        MINLEN:' $minLen
echo '        LEADING:' $qLead
echo '        TRAILING:' $qTrail
java -Xms4g -Xmx4g -jar /soft/bioinfo/linux_RH7/mugqic_pipelines-2.2.0/resources/software/trimmomatic/Trimmomatic-0.35/trimmomatic-0.35.jar PE \
             -threads $NB_THREADS -phred33 \
              ./$SAMPLE_NAME'_R1'.fastq.gz \
              ./$SAMPLE_NAME'_R2'.fastq.gz \
	      ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq.gz \
	      ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq.gz \
	      ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq.gz \
	      ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq.gz \
             ILLUMINACLIP:/soft/bioinfo/linux_RH7/mugqic_pipelines-2.2.0/resources/software/trimmomatic/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10:8:true \
             LEADING:$qLead TRAILING:$qTrail MINLEN:$minLen
echo ''

echo '    Number of R1 reads:'
echo '        Total: '`zcat ./$SAMPLE_NAME'_R1'.fastq.gz | wc -l`
echo '        Trimmed paired: '`zcat ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq.gz | wc -l`
echo '        Trimmed unpaired: '`zcat ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq.gz | wc -l`
# cat ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq > ./$SAMPLE_NAME'_R1'.trim.fastq
echo ''

echo '    Number of R2 reads:'
echo '        Total: '`zcat ./$SAMPLE_NAME'_R2'.fastq.gz | wc -l`
echo '        Trimmed paired: '`zcat ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq.gz | wc -l`
echo '        Trimmed unpaired: '`zcat ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq.gz | wc -l`
# cat ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq > ./$SAMPLE_NAME'_R2'.trim.fastq
echo ''

echo '    Deleting unused files...' `date "+%Y-%m-%d %H:%M:%S"`
rm -v ./$SAMPLE_NAME'_R1'.fastq.gz
rm -v ./$SAMPLE_NAME'_R2'.fastq.gz
echo ''

echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'
echo ''



if [[ $STRANDED = 'TRUE' ]]; then
    # Reverse Complement R1 and stranded jf database
    echo `date "+%Y-%m-%d %H:%M:%S":` 'Branch for stranded data...'
    echo '    Reverse complementing R1 fastq...'  `date "+%Y-%m-%d %H:%M:%S"`
    zcat ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq.gz | fastx_reverse_complement -z > ./$SAMPLE_NAME'_R1'.trimmed.paired.rc.fastq.gz
    zcat ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq.gz | fastx_reverse_complement -z > ./$SAMPLE_NAME'_R1'.trimmed.unpaired.rc.fastq.gz
    # python /u/laumontc/tsaPaper/fastq_reverse_complement.py ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq
    # python /u/laumontc/tsaPaper/fastq_reverse_complement.py ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq

    # echo '    Zipping fastq files...' `date "+%Y-%m-%d %H:%M:%S"`
    # gzip -v ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq
    # gzip -v ./$SAMPLE_NAME'_R1'.trimmed.paired.rc.fastq
    # gzip -v ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq
    # gzip -v ./$SAMPLE_NAME'_R1'.trimmed.unpaired.rc.fastq
    # gzip -v ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq
    # gzip -v ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq

    echo '    Copying trimmed fastq.gz files to outpath...' `date "+%Y-%m-%d %H:%M:%S"`
    cp -v ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq.gz $OUTPATH
    cp -v ./$SAMPLE_NAME'_R1'.trimmed.paired.rc.fastq.gz $OUTPATH
    cp -v ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq.gz $OUTPATH
    cp -v ./$SAMPLE_NAME'_R1'.trimmed.unpaired.rc.fastq.gz $OUTPATH
    cp -v ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq.gz $OUTPATH
    cp -v ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq.gz $OUTPATH

    echo '    Generating non-canonical jf database and other info...' `date "+%Y-%m-%d %H:%M:%S"`
    for l in "${LEN[@]}"
    do
	echo '        Length:' $l'...' `date "+%Y-%m-%d %H:%M:%S"`
	jellyfish count -m $l -s $HASH_SIZE -t $NB_THREADS -o $OUTPATH/$SAMPLE_NAME.trim.R1rc.$l.jf \
		  <(zcat ./$SAMPLE_NAME'_R1'.trimmed.paired.rc.fastq.gz) \
		  <(zcat ./$SAMPLE_NAME'_R1'.trimmed.unpaired.rc.fastq.gz) \
		  <(zcat ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq.gz) \
		  <(zcat ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq.gz)
	jellyfish histo $OUTPATH/$SAMPLE_NAME.trim.R1rc.$l.jf > $OUTPATH/$SAMPLE_NAME.trim.R1rc.$l.histo
	jellyfish info $OUTPATH/$SAMPLE_NAME.trim.R1rc.$l.jf > $OUTPATH/$SAMPLE_NAME.trim.R1rc.$l.info
    done

    # echo '    Deleting fastq.gz files...' `date "+%Y-%m-%d %H:%M:%S"`
    # rm -v ./*.fastq.gz
    # echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'

else
    # No reverse complement and unstranded jf database
    echo `date "+%Y-%m-%d %H:%M:%S":` 'Branch for unstranded data...'

    # echo '    Zipping fastq files...' `date "+%Y-%m-%d %H:%M:%S"`
    # gzip -v ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq
    # gzip -v ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq
    # gzip -v ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq
    # gzip -v ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq

    echo '    Copying trimmed fastq.gz files to outpath...' `date "+%Y-%m-%d %H:%M:%S"`
    cp -v ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq.gz $OUTPATH
    cp -v ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq.gz $OUTPATH
    cp -v ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq.gz $OUTPATH
    cp -v ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq.gz $OUTPATH

    echo '    Generating canonical jf database (-C) and other info...' `date "+%Y-%m-%d %H:%M:%S"`
    for l in "${LEN[@]}"
    do
	echo '        Length:' $l'...' `date "+%Y-%m-%d %H:%M:%S"`
	jellyfish count -m $l -s $HASH_SIZE -t $NB_THREADS -C -o $OUTPATH/$SAMPLE_NAME.trim.$l.jf \
		  <(zcat ./$SAMPLE_NAME'_R1'.trimmed.paired.fastq.gz) \
		  <(zcat ./$SAMPLE_NAME'_R1'.trimmed.unpaired.fastq.gz) \
		  <(zcat ./$SAMPLE_NAME'_R2'.trimmed.paired.fastq.gz) \
		  <(zcat ./$SAMPLE_NAME'_R2'.trimmed.unpaired.fastq.gz)
	jellyfish histo $OUTPATH/$SAMPLE_NAME.trim.$l.jf > $OUTPATH/$SAMPLE_NAME.trim.$l.histo
	jellyfish info $OUTPATH/$SAMPLE_NAME.trim.$l.jf > $OUTPATH/$SAMPLE_NAME.trim.$l.info
    done

    # echo '    Deleting fastq.gz files...' `date "+%Y-%m-%d %H:%M:%S"`
    # rm -v ./*.fastq.gz
    # echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'

fi

## Cleanup
echo 'Final cleanup' `date "+%Y-%m-%d %H:%M:%S"`
cd ..
rm -vrf $SAMPLE_NAME

echo `date "+%Y-%m-%d %H:%M:%S":` 'Done!'