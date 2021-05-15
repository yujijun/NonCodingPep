# Create back-up of .pyGeno file
# mkdir /u/laumontc/bck
# cp -r /u/laumontc/.pyGeno /u/laumontc/bck/.pyGeno

# Rename .pyGeno file
echo 'Copying .pyGeno to the RAM'
mv /u/laumontc/.pyGeno /u/laumontc/.pyGeno2
# Copy it to the RAM
mkdir -p /dev/shm/laumontc
cp -r /u/laumontc/.pyGeno2 /dev/shm/laumontc/.pyGeno
# Create symbolic link (pyGeno is expecting a .pyGeno file to run)
ln -s /dev/shm/laumontc/.pyGeno /u/laumontc/.pyGeno
echo 'Done!'

echo 'PT generation....'
python /u/laumontc/tsaPaper/scripts/cDNAperso.py -v GRCh38.88 -s plc07h103 -snp plc07h103mut -o /u/laumontc/tsaPaper/plc07h103/msRes/annotationFiles
echo 'Done!'

# Delete symbolic link
echo 'Restauring .pyGeno stuff...'
rm -f /u/laumontc/.pyGeno
# Rename your initial file as .pyGeno rather than .pyGeno2
mv /u/laumontc/.pyGeno2 /u/laumontc/.pyGeno
echo 'Done!'
