PATH_IN='/home/binsrv5/encode/'
PATH_OUT='/u/laumontc/tsaPaper/encode/natMed'

for t in $(\ls $PATH_IN)
do
    if [[ -d $PATH_IN/$t ]] ; then
	# echo $t
	echo $t >> $PATH_OUT/tissues.tmp
	echo $(ls $PATH_IN/$t/ | grep -E -o "SRR[0-9]+" | uniq | tr '\n' ' ') > $PATH_OUT/$t.srr
    fi
done

cat $PATH_OUT/*.srr > $PATH_OUT/tissueSRR.tmp
paste $PATH_OUT/tissues.tmp $PATH_OUT/tissueSRR.tmp > $PATH_OUT/tissues.txt

rm $PATH_OUT/*.srr $PATH_OUT/*.tmp

